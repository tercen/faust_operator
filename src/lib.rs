//! faust_operator — FAUST for Tercen, a drop-in for the R `faust_operator`.
//!
//! Rows are markers, columns are cells (optionally with the sample they belong to), y is the
//! value. The result is one row per cell: the FAUST phenotype label, one gate (0/1) per
//! selected marker, and the two coordinates of the annotation embedding.
//!
//! The algorithm is the `faust` crate (`crates/faust`): a port of FAUST 0.5.5 + scamp 0.5.0
//! checked stage by stage against the R packages. Samples are processed in parallel (the R
//! code loops over them); the search inside a sample is single-threaded, which is what makes
//! the result reproducible and identical to R's `threadNum = 1`.
pub mod context;
pub mod input;
pub mod output;
pub mod pagecache;
pub mod plots;
pub mod progress;
pub mod props;
pub mod tson;
pub mod upload;

use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result, anyhow};
use ndarray::Array2;
use tercen_rs::context::ContextBase;
use tercen_rs::{DevContext, TercenClient};

use faust::pipeline::{FaustSettings, Sample, run_faust};
use progress::Reporter;
use tson::TsonWriter;

const CHUNK: usize = 1_000_000;

pub fn init_tracing() {
    let filter = tracing_subscriber::EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info"));
    let _ = tracing_subscriber::fmt().with_env_filter(filter).try_init();
}

pub fn require_env(name: &str) -> Result<String> {
    std::env::var(name).map_err(|_| anyhow!("{name} is not set"))
}

pub async fn run(task_id: &str) -> Result<()> {
    tracing::info!("faust_operator starting (task_id={task_id})");
    let client = build_client().await?;
    let ctx = context::from_task_id(client, task_id).await?;
    execute(
        &ctx,
        Mode::Production {
            task_id: task_id.to_string(),
        },
    )
    .await
}

pub async fn run_dev(workflow_id: &str, step_id: &str) -> Result<()> {
    tracing::info!("faust_operator starting in dev mode ({workflow_id} / {step_id})");
    let client = build_client().await?;
    let ctx = DevContext::from_workflow_step(client, workflow_id, step_id)
        .await
        .map_err(|e| anyhow!("load workflow {workflow_id} / step {step_id}: {e}"))?;
    execute(
        &ctx,
        Mode::Dev {
            workflow_id: workflow_id.to_string(),
            step_id: step_id.to_string(),
        },
    )
    .await
}

enum Mode {
    Production {
        task_id: String,
    },
    Dev {
        workflow_id: String,
        step_id: String,
    },
}

async fn build_client() -> Result<Arc<TercenClient>> {
    let client = TercenClient::from_env()
        .await
        .map_err(|e| anyhow!("connect to Tercen: {e}"))?;
    tracing::info!("connected to Tercen");
    Ok(Arc::new(client))
}

async fn execute(ctx: &ContextBase, mode: Mode) -> Result<()> {
    let t_start = Instant::now();
    tracing::info!(
        workflow = ctx.workflow_id(),
        step = ctx.step_id(),
        namespace = ctx.namespace(),
        "context loaded"
    );
    let rep = match &mode {
        Mode::Production { task_id } => Reporter::spawn(Arc::clone(ctx.client()), task_id.clone()),
        Mode::Dev { .. } => Reporter::silent(),
    };
    let s = props::read(ctx)?;

    let n_values = input::cell_count(ctx).await?;
    let markers = input::row_labels(ctx).await?;
    let p = markers.len();
    if p == 0 {
        anyhow::bail!("the projection has no rows: put the markers on rows");
    }
    let n_cells = n_values / p;
    tracing::info!(n_values, n_cells, markers = p, "projection");

    // The sample of every cell: the `filename` column factor as the R operator, or the
    // factor named by `sample_factor`. One sample when neither is projected.
    let groups = sample_groups(ctx, &s, n_cells, &rep).await?;
    let n_samples = groups.iter().copied().max().map_or(1, |m| m + 1);

    rep.at(0, "Reading the crosstab");
    let data = gather(ctx, n_values, n_cells, p, &rep).await?;

    // Split into per-sample, column-major matrices, remembering the cells' .ci.
    let mut cells_of: Vec<Vec<usize>> = vec![Vec::new(); n_samples];
    for (ci, &g) in groups.iter().enumerate() {
        cells_of[g].push(ci);
    }
    let samples: Vec<Sample> = cells_of
        .iter()
        .enumerate()
        .map(|(g, cells)| Sample {
            name: format!("sample_{g}"),
            cols: (0..p)
                .map(|j| cells.iter().map(|&ci| data[ci * p + j]).collect())
                .collect(),
        })
        .collect();
    drop(data);
    for (g, smp) in samples.iter().enumerate() {
        if smp.cols[0].len() < 25 {
            anyhow::bail!(
                "sample {g} has {} cells; FAUST needs at least 25 per sample (minimumClusterSize)",
                smp.cols[0].len()
            );
        }
    }

    rep.at(progress::FIT.0, format!("FAUST: {n_samples} samples x {p} markers"));
    let settings = FaustSettings {
        selection_quantile: s.selection_quantile,
        depth_score_threshold: s.depth_score_threshold,
        seed: s.seed,
        threads: s.threads,
    };
    let t = Instant::now();
    let rep2 = rep.clone();
    let run = run_faust(&samples, &markers, &settings, &mut |msg: &str| {
        rep2.info(msg.to_string());
    })?;
    tracing::info!(
        secs = format!("{:.1}", t.elapsed().as_secs_f64()),
        selected = run.selected.len(),
        "faust done"
    );
    let selected_names: Vec<String> = run.selected.iter().map(|&j| markers[j].clone()).collect();
    rep.info(format!(
        "FAUST selected {} of {p} markers: {}",
        run.selected.len(),
        selected_names.join(", ")
    ));

    // ---- annotation embedding
    rep.at(progress::band(progress::FIT, 3, 4), "Annotation embedding (UMAP)");
    let embed_samples: Vec<faust::embed::EmbedSample> = (0..n_samples)
        .map(|g| faust::embed::EmbedSample {
            exprs: &run.selected_exprs[g],
            ann: &run.ann[g],
        })
        .collect();
    let prepared = faust::embed::prepare(&embed_samples, &run.levels, &run.score_mat, s.seed as u64);
    let em = &prepared.embedding_markers;
    let (x, y) = if em.is_empty() || n_cells <= 15 {
        rep.info("too few markers or cells for the embedding; umap_1/umap_2 are NaN".to_string());
        (vec![f64::NAN; n_cells], vec![f64::NAN; n_cells])
    } else {
        let flat: Vec<f64> = prepared
            .rows
            .iter()
            .flat_map(|r| em.iter().map(move |&j| r[j]))
            .collect();
        let m = Array2::from_shape_vec((n_cells, em.len()), flat)
            .map_err(|e| anyhow!("shape the embedding matrix: {e}"))?;
        let umap = umaprs::UMAP::new()
            .n_neighbors(15)
            .min_dist(0.2)
            .random_state(s.seed as u64)
            .threads(s.threads);
        let e = umap.fit_transform(&m);
        // rows are in sample order (cells_of concatenated); map back to .ci
        let mut x = vec![f64::NAN; n_cells];
        let mut y = vec![f64::NAN; n_cells];
        let mut k = 0;
        for cells in &cells_of {
            for &ci in cells {
                x[ci] = e[[k, 0]];
                y[ci] = e[[k, 1]];
                k += 1;
            }
        }
        (x, y)
    };

    // ---- per-cell columns in .ci order
    let mut labels = vec![String::new(); n_cells];
    let mut gates: Vec<(String, Vec<f64>)> = selected_names
        .iter()
        .map(|n| (n.clone(), vec![f64::NAN; n_cells]))
        .collect();
    for (g, cells) in cells_of.iter().enumerate() {
        for (k, &ci) in cells.iter().enumerate() {
            labels[ci] = run.labels[g][k].clone();
            for (j, lvl) in run.ann[g][k].iter().enumerate() {
                // The R operator maps "-" to 0 and "+" to 1; other levels (Dim, Bright…) are NA.
                gates[j].1[ci] = if run.levels[j] == 2 {
                    match lvl {
                        1 => 0.0,
                        2 => 1.0,
                        _ => f64::NAN,
                    }
                } else {
                    f64::NAN
                };
            }
        }
    }

    // ---- histograms
    let pngs = if s.histograms {
        plots::marker_histograms(&run, &selected_names, s.seed as u64)?
    } else {
        Vec::new()
    };

    rep.at(progress::WRITE.0, "Writing the result");
    let work_root =
        std::env::temp_dir().join(format!("faust_op_{}_{}", ctx.workflow_id(), ctx.step_id()));
    std::fs::create_dir_all(&work_root)
        .with_context(|| format!("create {}", work_root.display()))?;
    let _guard = TempDirGuard(work_root.clone());
    let result_path = work_root.join("result.tson");
    {
        let f = std::fs::File::create(&result_path)
            .with_context(|| format!("create {}", result_path.display()))?;
        let w = std::io::BufWriter::with_capacity(4 << 20, pagecache::Releasing::new(f, 256 << 20));
        let mut w = TsonWriter::new(w)?;
        output::write_result(
            &mut w,
            &table_name(ctx),
            &output::CellTable {
                namespace: ctx.namespace(),
                umap_1: &x,
                umap_2: &y,
                labels: &labels,
                gates: &gates,
            },
            &pngs,
        )?;
    }
    let bytes = std::fs::metadata(&result_path)?.len();
    tracing::info!(bytes, "result written");
    pagecache::release_path(&result_path);

    rep.at(progress::UPLOAD.0, "Uploading the result");
    match mode {
        Mode::Production { task_id } => {
            upload::save_production(ctx, &task_id, &result_path, &rep).await?
        }
        Mode::Dev {
            workflow_id,
            step_id,
        } => {
            let saved = upload::save_dev(ctx, &workflow_id, &step_id, &result_path).await?;
            tracing::info!(task_id = saved.task_id, "dev result saved");
        }
    }
    rep.at(100, "Done");
    tracing::info!(
        total_secs = format!("{:.1}", t_start.elapsed().as_secs_f64()),
        peak_rss_kb = peak_rss_kb().unwrap_or(0),
        "done"
    );
    Ok(())
}

async fn sample_groups(
    ctx: &ContextBase,
    s: &props::Settings,
    n_cells: usize,
    rep: &Reporter,
) -> Result<Vec<usize>> {
    let names = ctx
        .cnames()
        .await
        .map_err(|e| anyhow!("read column factor names: {e}"))?;
    let wanted = if s.sample_factor.is_empty() {
        match names
            .iter()
            .find(|n| n.as_str() == "filename" || n.ends_with(".filename"))
        {
            Some(n) => n.clone(),
            None => {
                rep.info(format!(
                    "no 'filename' column factor is projected (columns are {names:?}); all \
                     {n_cells} cells are one sample. Set sample_factor to split by sample."
                ));
                return Ok(vec![0; n_cells]);
            }
        }
    } else {
        s.sample_factor.clone()
    };
    let g = input::column_groups(ctx, &wanted).await?;
    if g.len() != n_cells {
        anyhow::bail!(
            "the column factor '{wanted}' has {} values for {n_cells} cells",
            g.len()
        );
    }
    let n_groups = g.iter().copied().max().map_or(0, |m| m + 1);
    if n_groups == n_cells && n_cells > 1 {
        rep.info(format!(
            "the sample factor '{wanted}' has one value per cell; all cells count as one sample"
        ));
        return Ok(vec![0; n_cells]);
    }
    Ok(g)
}

/// Gather the crosstab into a cell-by-marker matrix (cell-major), scattered by `(.ri, .ci)`.
async fn gather(
    ctx: &ContextBase,
    n_values: usize,
    n_cells: usize,
    p: usize,
    rep: &Reporter,
) -> Result<Vec<f64>> {
    let mut data = vec![f64::NAN; n_cells * p];
    let mut seen = 0usize;
    let mut out_of_range = 0usize;
    input::for_each_chunk(ctx, &[".ri", ".ci", ".y"], n_values, CHUNK, |c| {
        for k in 0..c.y.len() {
            let (ri, ci) = (c.ri[k] as usize, c.ci[k] as usize);
            if ri >= p || ci >= n_cells {
                out_of_range += 1;
                continue;
            }
            data[ci * p + ri] = c.y[k];
        }
        seen += c.len();
        rep.at(
            progress::band(progress::READ, seen, n_values.max(1)),
            format!("Read {seen} of {n_values}"),
        );
        Ok(())
    })
    .await?;
    if out_of_range > 0 {
        anyhow::bail!(
            "{out_of_range} values fell outside the {n_cells} x {p} crosstab the schema \
             described; the projection changed under the run"
        );
    }
    let missing = data.iter().filter(|v| !v.is_finite()).count();
    if missing > 0 {
        anyhow::bail!(
            "{missing} of {} cell-marker pairs are missing or not finite. FAUST needs a complete \
             matrix: every cell must have every marker.",
            data.len()
        );
    }
    Ok(data)
}

struct TempDirGuard(std::path::PathBuf);
impl Drop for TempDirGuard {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn table_name(ctx: &ContextBase) -> String {
    format!("{}_{}", ctx.step_id(), ctx.qt_hash())
}

fn peak_rss_kb() -> Option<u64> {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1)?.parse().ok())
        })
}

//! Operator properties. The three the R operator has (`selectionQuantile`,
//! `depthScoreThreshold`, `seed`) with the same names and defaults, plus `sample_factor` (the
//! R operator hard-codes `filename`) and `threads`.
use anyhow::{Result, bail};
use tercen_rs::PropertyReader;
use tercen_rs::context::ContextBase;

#[derive(Debug, Clone)]
pub struct Settings {
    /// `selectionQuantile`: a marker is kept when this quantile of its depth score across
    /// samples is at least `depthScoreThreshold`.
    pub selection_quantile: f64,
    /// `depthScoreThreshold`.
    pub depth_score_threshold: f64,
    /// `seed`: `seedValue` of FAUST — the forest seed is derived from it the way R does.
    pub seed: u32,
    /// `sample_factor`: the column factor that names the sample. Empty means `filename` if
    /// projected, else one sample.
    pub sample_factor: String,
    /// `threads`: samples are processed in parallel; 0 uses every core the task has.
    pub threads: usize,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            selection_quantile: 0.5,
            depth_score_threshold: 0.01,
            seed: 123,
            sample_factor: String::new(),
            threads: 0,
        }
    }
}

pub fn read(ctx: &ContextBase) -> Result<Settings> {
    let pr = PropertyReader::from_operator_settings(ctx.operator_settings());
    let d = Settings::default();
    let num = |name: &str, dflt: f64| -> Result<f64> {
        let raw = pr.get_string(name, "").trim().to_string();
        if raw.is_empty() || raw.eq_ignore_ascii_case("null") {
            return Ok(dflt);
        }
        raw.parse::<f64>()
            .map_err(|_| anyhow::anyhow!("property '{name}' should be a number, got '{raw}'"))
    };
    let text = |name: &str| pr.get_string(name, "").trim().to_string();
    let selection_quantile = num("selectionQuantile", d.selection_quantile)?;
    if !(0.0..=1.0).contains(&selection_quantile) {
        bail!("property 'selectionQuantile' must be between 0 and 1, got {selection_quantile}");
    }
    let depth_score_threshold = num("depthScoreThreshold", d.depth_score_threshold)?;
    if !(0.0..=1.0).contains(&depth_score_threshold) {
        bail!(
            "property 'depthScoreThreshold' must be between 0 and 1, got {depth_score_threshold}"
        );
    }
    let seed = num("seed", d.seed as f64)?;
    if !(1.0..=4_294_967_295.0).contains(&seed) {
        bail!("property 'seed' must be a positive integer below 2^32, got {seed}");
    }
    let s = Settings {
        selection_quantile,
        depth_score_threshold,
        seed: seed as u32,
        sample_factor: text("sample_factor"),
        threads: num("threads", 0.0)?.max(0.0) as usize,
    };
    tracing::info!(?s, "properties");
    Ok(s)
}

#[cfg(test)]
mod tests {
    /// Tercen books `intercept × n_cols^1 × n_rows^0.2 + 1.5 × offset` MB: cells (the column
    /// table) times a mild growth with markers (the row table). Peak RAM follows cells — the UMAP
    /// of every cell, the gathered matrix and the per-cell output — far more than values: 8 → 23
    /// markers on the same 1.42 M cells raised it 19 %. Measured `stats_d_actual_ram_peak` (MB):
    /// tercen.com 1.1.2 at 6 markers (60k, 75k, 200k cells), Studio 1.2.0 at 8, 23 and 31 markers.
    /// Peaks are in 10^6 bytes and Tercen books MiB, so the margins here are ~5 % conservative.
    /// The previous per-value model booked 3.6× the 23-marker run's peak.
    #[test]
    fn memory_model_books_measured_peaks_by_cells() {
        let m: serde_json::Value =
            serde_json::from_str(include_str!("../memory_model.json")).unwrap();
        let f = m["features"].as_array().unwrap();
        assert_eq!(
            (f[0]["name"].as_str(), f[1]["name"].as_str()),
            (Some("n_cols"), Some("n_rows"))
        );
        let book = |cells: f64, markers: f64| {
            m["intercept"].as_f64().unwrap()
                * cells.powf(f[0]["exponent"].as_f64().unwrap())
                * markers.powf(f[1]["exponent"].as_f64().unwrap())
                + 1.5 * m["offset"].as_f64().unwrap()
        };
        for (cells, markers, peak) in [
            (60_000.0, 6.0, 170.0),
            (75_000.0, 6.0, 212.0),
            (200_000.0, 6.0, 518.0),
            (1_420_000.0, 8.0, 3_944.0),
            (1_420_000.0, 23.0, 4_590.0),
            (1_420_000.0, 31.0, 4_880.0),
        ] {
            let b = book(cells, markers);
            assert!(
                b >= 1.25 * peak,
                "{cells} x {markers}: books {b:.0} MB for a {peak} MB peak"
            );
            if cells >= 1e6 {
                assert!(
                    b <= 2.0 * peak,
                    "{cells} x {markers}: books {b:.0} MB, over twice the peak"
                );
            }
        }
    }

    use super::*;
    use serde_json::Value;

    #[test]
    fn manifest_declares_what_the_code_reads() {
        let m: Value = serde_json::from_str(include_str!("../operator.json")).unwrap();
        let props: Vec<&Value> = m["properties"].as_array().unwrap().iter().collect();
        let by_name = |n: &str| {
            props
                .iter()
                .find(|p| p["name"] == n)
                .unwrap_or_else(|| panic!("property '{n}' is not declared in operator.json"))
        };
        let d = Settings::default();
        assert_eq!(
            by_name("selectionQuantile")["defaultValue"],
            d.selection_quantile
        );
        assert_eq!(
            by_name("depthScoreThreshold")["defaultValue"],
            d.depth_score_threshold
        );
        assert_eq!(by_name("seed")["defaultValue"], d.seed as f64);
        assert_eq!(by_name("sample_factor")["defaultValue"], "");
        assert_eq!(by_name("threads")["defaultValue"], 0.0);
        let allowed = [
            "OperatorSpec",
            "CrosstabSpec",
            "MetaFactor",
            "AxisSpec",
            "OperatorJoinSpec",
            "JoinOperator",
            "ColumnPair",
            "TableRelation",
            "Attribute",
            "Pair",
        ];
        fn walk(v: &Value, allowed: &[&str]) {
            match v {
                Value::Object(o) => {
                    if let Some(Value::String(k)) = o.get("kind") {
                        assert!(allowed.contains(&k.as_str()), "invented kind '{k}'");
                    }
                    o.values().for_each(|c| walk(c, allowed));
                }
                Value::Array(a) => a.iter().for_each(|c| walk(c, allowed)),
                _ => {}
            }
        }
        walk(&m["operatorSpec"], &allowed);
        let c = m["container"].as_str().unwrap();
        let tag = c.rsplit_once(':').map(|(_, t)| t).unwrap_or("");
        assert!(
            tag.split('.').count() == 3
                && tag
                    .split('.')
                    .all(|p| !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit())),
            "container must pin an exact x.y.z version tag"
        );
        assert_eq!(m["communicationProtocol"], "grpc");
    }
}

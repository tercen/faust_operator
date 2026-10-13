//! The result tables, written as an `OperatorResult` TSON stream.
//!
//! Table 1, one row per cell (keyed by `.ci`): `umap_1`, `umap_2`, `label`, one `<marker>_gate`
//! per selected marker — the R operator's columns, with the same types (the gates are doubles
//! there and stay doubles so a workflow built on 0.x keeps working).
//!
//! Table 2 (when there are histograms), one row per selected marker keyed by `.ri`: `filename`,
//! `mimetype`, `checksum`, `.content` (base64 PNG) — the shape `tim::plot_file_to_df` gives the
//! R operator, so the UI shows them under "Files". The R operator joined its PNG table on
//! nothing (a cross join); keyed on `.ri` the image sits on the marker's row instead.
//!
//! `joinOperators` stays empty: that is the branch where the server joins every table itself
//! by its `.ci`/`.ri` columns. Declaring joins switches it to the explicit branch, where the
//! per-cell table was no longer attached (1.0.0 dropped the PNGs for that reason).
use std::io::Write;

use anyhow::{Result, anyhow};

use crate::tson::TsonWriter;

pub struct ColSpec<'a> {
    pub name: &'a str,
    pub ty: &'a str, // "double" | "int32" | "string"
}

fn table_header<W: Write>(
    w: &mut TsonWriter<W>,
    table_name: &str,
    n_rows: usize,
    n_cols: usize,
) -> Result<()> {
    w.map(4)?;
    w.key("kind")?;
    w.str("Table")?;
    w.key("nRows")?;
    w.i32(i32::try_from(n_rows).map_err(|_| {
        anyhow!(
            "the result would have {n_rows} rows, more than a Tercen table can hold (i32::MAX). \
             Project fewer cells, or split the step."
        )
    })?)?;
    w.key("properties")?;
    w.map(4)?;
    w.key("kind")?;
    w.str("TableProperties")?;
    w.key("name")?;
    w.str(table_name)?;
    w.key("sortOrder")?;
    w.list(0)?;
    w.key("ascending")?;
    w.bool(false)?;
    w.key("columns")?;
    w.list(n_cols)?;
    Ok(())
}

fn column_header<W: Write>(w: &mut TsonWriter<W>, c: &ColSpec, n_rows: usize) -> Result<()> {
    w.map(6)?;
    w.key("kind")?;
    w.str("Column")?;
    w.key("name")?;
    w.str(c.name)?;
    w.key("type")?;
    w.str(c.ty)?;
    w.key("nRows")?;
    w.i32(n_rows as i32)?;
    w.key("size")?;
    w.i32(n_rows as i32)?;
    w.key("values")?;
    Ok(())
}

pub struct CellTable<'a> {
    pub namespace: &'a str,
    pub umap_1: &'a [f64],
    pub umap_2: &'a [f64],
    pub labels: &'a [String],
    /// `(marker, 0/1 per cell)` in the order the columns are written.
    pub gates: &'a [(String, Vec<f64>)],
}

/// The histogram table: parallel vectors, one entry per selected marker.
pub struct HistTable {
    /// Row index of the marker in the crosstab.
    pub ri: Vec<i32>,
    pub filenames: Vec<String>,
    /// md5 of the PNG bytes, hex.
    pub checksums: Vec<String>,
    /// The PNG, base64.
    pub contents: Vec<String>,
}

pub fn write_result<W: Write>(
    w: &mut TsonWriter<W>,
    table_name: &str,
    cells: &CellTable,
    hists: Option<&HistTable>,
) -> Result<()> {
    let n = cells.umap_1.len();
    let ns = cells.namespace;
    let hists = hists.filter(|h| !h.ri.is_empty());
    w.map(3)?;
    w.key("kind")?;
    w.str("OperatorResult")?;
    w.key("tables")?;
    w.list(if hists.is_some() { 2 } else { 1 })?;

    // ---- table 1: cells
    let names: Vec<String> = [
        format!("{ns}.umap_1"),
        format!("{ns}.umap_2"),
        format!("{ns}.label"),
    ]
    .into_iter()
    .chain(cells.gates.iter().map(|(m, _)| format!("{ns}.{m}_gate")))
    .chain(std::iter::once(".ci".to_string()))
    .collect();
    table_header(w, table_name, n, names.len())?;
    let mut k = 0;
    column_header(
        w,
        &ColSpec {
            name: &names[k],
            ty: "double",
        },
        n,
    )?;
    w.f64_list(cells.umap_1)?;
    k += 1;
    column_header(
        w,
        &ColSpec {
            name: &names[k],
            ty: "double",
        },
        n,
    )?;
    w.f64_list(cells.umap_2)?;
    k += 1;
    column_header(
        w,
        &ColSpec {
            name: &names[k],
            ty: "string",
        },
        n,
    )?;
    w.str_list(cells.labels)?;
    k += 1;
    for (_, g) in cells.gates {
        column_header(
            w,
            &ColSpec {
                name: &names[k],
                ty: "double",
            },
            n,
        )?;
        w.f64_list(g)?;
        k += 1;
    }
    column_header(
        w,
        &ColSpec {
            name: &names[k],
            ty: "int32",
        },
        n,
    )?;
    w.i32_list(&(0..n as i32).collect::<Vec<_>>())?;

    // ---- table 2: histograms, keyed by .ri
    if let Some(h) = hists {
        let m = h.ri.len();
        table_header(w, &format!("{table_name}_histograms"), m, 5)?;
        let mimetypes = vec!["image/png".to_string(); m];
        for (name, values) in [
            (format!("{ns}.filename"), &h.filenames),
            (format!("{ns}.mimetype"), &mimetypes),
            (format!("{ns}.checksum"), &h.checksums),
            (".content".to_string(), &h.contents),
        ] {
            column_header(
                w,
                &ColSpec {
                    name: &name,
                    ty: "string",
                },
                m,
            )?;
            w.str_list(values)?;
        }
        column_header(
            w,
            &ColSpec {
                name: ".ri",
                ty: "int32",
            },
            m,
        )?;
        w.i32_list(&h.ri)?;
    }

    // No joinOperators: the server joins each table by its .ci / .ri column.
    w.key("joinOperators")?;
    w.list(0)?;
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_cell_table_decodes_with_the_expected_columns() {
        let mut w = TsonWriter::new(Vec::new()).unwrap();
        let gates = vec![("CD3".to_string(), vec![1.0, 0.0])];
        write_result(
            &mut w,
            "t",
            &CellTable {
                namespace: "ds0",
                umap_1: &[1.0, 2.0],
                umap_2: &[3.0, 4.0],
                labels: &["CD3+".to_string(), "CD3-".to_string()],
                gates: &gates,
            },
            None,
        )
        .unwrap();
        let bytes = w.into_inner();
        let v = rustson::decode(std::io::Cursor::new(&bytes[..])).unwrap();
        let rustson::Value::MAP(m) = v else { panic!() };
        let rustson::Value::LST(tables) = &m["tables"] else {
            panic!()
        };
        let rustson::Value::MAP(t) = &tables[0] else {
            panic!()
        };
        let rustson::Value::LST(cols) = &t["columns"] else {
            panic!()
        };
        let names: Vec<String> = cols
            .iter()
            .map(|c| {
                let rustson::Value::MAP(c) = c else { panic!() };
                let rustson::Value::STR(n) = &c["name"] else {
                    panic!()
                };
                n.clone()
            })
            .collect();
        assert_eq!(
            names,
            [
                "ds0.umap_1",
                "ds0.umap_2",
                "ds0.label",
                "ds0.CD3_gate",
                ".ci"
            ]
        );
    }

    #[test]
    fn the_histogram_table_is_keyed_by_ri() {
        let mut w = TsonWriter::new(Vec::new()).unwrap();
        let gates = vec![("CD3".to_string(), vec![1.0, 0.0])];
        let h = HistTable {
            ri: vec![2],
            filenames: vec!["hist_CD3.png".into()],
            checksums: vec!["abc".into()],
            contents: vec!["iVBORw0KGgo=".into()],
        };
        write_result(
            &mut w,
            "t",
            &CellTable {
                namespace: "ds0",
                umap_1: &[1.0, 2.0],
                umap_2: &[3.0, 4.0],
                labels: &["CD3+".to_string(), "CD3-".to_string()],
                gates: &gates,
            },
            Some(&h),
        )
        .unwrap();
        let v = rustson::decode(std::io::Cursor::new(&w.into_inner()[..])).unwrap();
        let rustson::Value::MAP(m) = v else { panic!() };
        let rustson::Value::LST(tables) = &m["tables"] else {
            panic!()
        };
        assert_eq!(tables.len(), 2);
        let rustson::Value::LST(joins) = &m["joinOperators"] else {
            panic!()
        };
        assert!(joins.is_empty());
        let rustson::Value::MAP(t) = &tables[1] else {
            panic!()
        };
        let rustson::Value::LST(cols) = &t["columns"] else {
            panic!()
        };
        let names: Vec<String> = cols
            .iter()
            .map(|c| {
                let rustson::Value::MAP(c) = c else { panic!() };
                let rustson::Value::STR(n) = &c["name"] else {
                    panic!()
                };
                n.clone()
            })
            .collect();
        assert_eq!(
            names,
            [
                "ds0.filename",
                "ds0.mimetype",
                "ds0.checksum",
                ".content",
                ".ri"
            ]
        );
    }
}

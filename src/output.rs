//! The result tables, written as an `OperatorResult` TSON stream.
//!
//! Table 1, one row per cell (keyed by `.ci`): `umap_1`, `umap_2`, `label`, one `<marker>_gate`
//! per selected marker — the R operator's columns, with the same types (the gates are doubles
//! there and stay doubles so a workflow built on 0.x keeps working).
//!
//! The R operator's third table (the threshold histograms as PNGs, joined on nothing) is not
//! produced: an empty-key join is a cross join against the crosstab, and the platform dropped
//! the main table when the result declared one.
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

/// The whole result: the per-cell table and, when `pngs` is non-empty, the histogram table.
pub struct CellTable<'a> {
    pub namespace: &'a str,
    pub umap_1: &'a [f64],
    pub umap_2: &'a [f64],
    pub labels: &'a [String],
    /// `(marker, 0/1 per cell)` in the order the columns are written.
    pub gates: &'a [(String, Vec<f64>)],
}

pub fn write_result<W: Write>(
    w: &mut TsonWriter<W>,
    table_name: &str,
    cells: &CellTable,
) -> Result<()> {
    let n = cells.umap_1.len();
    let ns = cells.namespace;
    w.map(3)?;
    w.key("kind")?;
    w.str("OperatorResult")?;
    w.key("tables")?;
    w.list(1)?;

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

    // One relation only, keyed by .ci: no joinOperators (a second table joined on nothing is
    // a cross join against the crosstab, and the platform did not persist the main table when
    // one was declared).
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
}

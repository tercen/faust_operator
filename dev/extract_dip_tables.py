#!/usr/bin/env python3
"""Extract FAUST's dip null tables from the reference C++ into little-endian f64 files.

Run from anywhere:  python3 -I dev/extract_dip_tables.py

Reads (relative to the repo root, resolved from this script's location)
  crates/faust/reference_cpp_faust/singleDip.cpp   -> _dip_0 .. _dip_20 (21 x 26, diptest qDiptab),
                                                      _dip1_row_names_ (21 sample sizes),
                                                      _dip1_col_names_ (26 probabilities)
  crates/faust/reference_cpp_faust/doubleDip.cpp   -> _dip2v_0_ .. _dip2v_499_ (500 x 100)
  crates/faust/reference_cpp_faust/tripleDip.cpp   -> _dip3v_0_ .. _dip3v_399_ (400 x 100)

and writes into crates/faust/src/dip_tables/ (embedded by `include_bytes!` in src/dip.rs):
  dip1_crit.f64le   21*26  f64, row-major, rows in the order of `_dip1_crit_vals_`
  dip1_rows.f64le   21     f64 (the sample sizes, exactly representable integers)
  dip1_cols.f64le   26     f64 (the probabilities)
  dip2_crit.f64le   500*100 f64, row-major, rows in the order of `_dip2_critical_values_`
                    (row N-1 is for sample size N)
  dip3_crit.f64le   400*100 f64, row-major, rows in the order of `_dip3_critical_values_`

Every number is parsed with Python's correctly rounded `float()`, which yields the same
IEEE double a conforming C++ compiler produces for the same decimal literal.
"""

import re
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "crates" / "faust" / "reference_cpp_faust"
OUT = ROOT / "crates" / "faust" / "src" / "dip_tables"

VEC_RE = re.compile(
    r"std::vector<(?:double|long|int)>\s+(_\w+)\s*=\s*\{([^{}]*)\}\s*;", re.S
)
AGG_RE = re.compile(
    r"std::vector<std::vector<double>>\s+(_\w+)\s*=\s*\{([^{}]*)\}\s*;", re.S
)


def parse(path):
    text = path.read_text()
    vecs = {}
    for name, body in VEC_RE.findall(text):
        items = [t.strip() for t in body.split(",") if t.strip()]
        vecs[name] = [float(t) for t in items]
    aggs = {}
    for name, body in AGG_RE.findall(text):
        aggs[name] = [t.strip() for t in body.split(",") if t.strip()]
    return vecs, aggs


def write_f64(path, values):
    path.write_bytes(struct.pack("<%dd" % len(values), *values))
    print("%-18s %7d f64" % (path.name, len(values)))


def matrix(vecs, aggs, agg_name, nrow, ncol):
    rows = aggs[agg_name]
    if len(rows) != nrow:
        sys.exit("%s: expected %d rows, found %d" % (agg_name, nrow, len(rows)))
    flat = []
    for r in rows:
        v = vecs[r]
        if len(v) != ncol:
            sys.exit("%s: expected %d values, found %d" % (r, ncol, len(v)))
        flat.extend(v)
    return flat


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    vecs, aggs = parse(REF / "singleDip.cpp")
    write_f64(OUT / "dip1_crit.f64le", matrix(vecs, aggs, "_dip1_crit_vals_", 21, 26))
    rows = vecs["_dip1_row_names_"]
    cols = vecs["_dip1_col_names_"]
    if len(rows) != 21 or len(cols) != 26:
        sys.exit("singleDip.cpp: unexpected row/col name lengths")
    if any(r != int(r) for r in rows):
        sys.exit("singleDip.cpp: non-integer sample size")
    write_f64(OUT / "dip1_rows.f64le", rows)
    write_f64(OUT / "dip1_cols.f64le", cols)

    vecs, aggs = parse(REF / "doubleDip.cpp")
    write_f64(OUT / "dip2_crit.f64le", matrix(vecs, aggs, "_dip2_critical_values_", 500, 100))

    vecs, aggs = parse(REF / "tripleDip.cpp")
    write_f64(OUT / "dip3_crit.f64le", matrix(vecs, aggs, "_dip3_critical_values_", 400, 100))


if __name__ == "__main__":
    main()

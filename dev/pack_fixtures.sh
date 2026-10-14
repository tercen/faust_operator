#!/bin/sh
# Rebuild crates/faust/tests/parity/fixtures.tar.gz from the fixtures/ directory
# (after dev/make_fixtures.R & co). Reproducible: sorted names, fixed owner and mtime, gzip -n.
set -eu
cd "$(dirname "$0")/../crates/faust/tests/parity"
tar --sort=name --owner=0 --group=0 --numeric-owner --mtime='2026-01-01 00:00Z' \
    --exclude='.unpacked' -cf - fixtures | gzip -9 -n > fixtures.tar.gz
ls -l fixtures.tar.gz

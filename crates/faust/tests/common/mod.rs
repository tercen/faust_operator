//! Where the parity fixtures are.
//!
//! They are tracked as one file, `tests/parity/fixtures.tar.gz` (820 files, 17 MB unpacked):
//! Tercen's library install imports every file of an operator repo as a document, and one
//! archive costs it one document instead of 820. `fixtures()` unpacks the archive next to
//! itself on first use and returns the directory; the unpacked tree is git-ignored. Rebuild the
//! archive with `dev/pack_fixtures.sh` after regenerating fixtures.
#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::sync::OnceLock;

pub fn fixtures() -> PathBuf {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    DIR.get_or_init(unpack).clone()
}

fn unpack() -> PathBuf {
    let parity = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/parity");
    let dir = parity.join("fixtures");
    let archive = parity.join("fixtures.tar.gz");
    if is_fresh(&dir, &archive) {
        return dir;
    }
    let tmp = parity.join(format!("fixtures.tmp-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).unwrap();
    let file = std::fs::File::open(&archive)
        .unwrap_or_else(|e| panic!("{}: {e} (run dev/pack_fixtures.sh)", archive.display()));
    tar::Archive::new(flate2::read::GzDecoder::new(file))
        .unpack(&tmp)
        .unwrap_or_else(|e| panic!("unpacking {}: {e}", archive.display()));
    std::fs::File::create(tmp.join("fixtures").join(".unpacked")).unwrap();
    if !is_fresh(&dir, &archive) {
        let _ = std::fs::remove_dir_all(&dir);
    }
    // Atomic publish; a concurrent test binary may have published first, which is fine.
    match std::fs::rename(tmp.join("fixtures"), &dir) {
        Ok(()) => {}
        Err(_) if is_fresh(&dir, &archive) => {}
        Err(e) => panic!("{}: {e}", dir.display()),
    }
    let _ = std::fs::remove_dir_all(&tmp);
    dir
}

/// The unpacked tree exists and is at least as new as the archive.
fn is_fresh(dir: &Path, archive: &Path) -> bool {
    let modified = |p: &Path| std::fs::metadata(p).and_then(|m| m.modified()).ok();
    match (modified(&dir.join(".unpacked")), modified(archive)) {
        (Some(unpacked), Some(packed)) => unpacked >= packed,
        _ => false,
    }
}

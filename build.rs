use std::env;
use std::fs;
use std::path::PathBuf;

fn main() {
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let proj_db_dest = out_dir.join("proj.db");

    // Only copy if not already present (incremental builds).
    if proj_db_dest.exists() {
        println!("cargo:rerun-if-changed=build.rs");
        return;
    }

    // Find proj.db: check PROJ_DB_PATH env var, then standard system locations.
    let candidates = [
        env::var("PROJ_DB_PATH").ok(),
        Some("/usr/share/proj/proj.db".to_string()),
        Some("/usr/local/share/proj/proj.db".to_string()),
        Some("/opt/homebrew/share/proj/proj.db".to_string()),
    ];

    let proj_db_src = candidates
        .iter()
        .filter_map(|c| c.as_ref())
        .map(PathBuf::from)
        .find(|p| p.exists());

    match proj_db_src {
        Some(src) => {
            fs::copy(&src, &proj_db_dest).unwrap_or_else(|e| {
                panic!("Failed to copy proj.db from {} to {}: {}", src.display(), proj_db_dest.display(), e)
            });
            eprintln!("Embedded proj.db from {}", src.display());
        }
        None => {
            panic!(
                "proj.db not found. Install proj-data (apt-get install proj-data) \
                 or set PROJ_DB_PATH to the path of proj.db."
            );
        }
    }

    println!("cargo:rerun-if-changed=build.rs");
    if let Ok(path) = env::var("PROJ_DB_PATH") {
        println!("cargo:rerun-if-changed={}", path);
    }
}

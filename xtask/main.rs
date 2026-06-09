use std::{
    env,
    error::Error,
    fs,
    path::{Path, PathBuf},
    process::{exit, Command},
};

fn main() {
    if let Err(err) = run() {
        eprintln!("xtask: {err}");
        exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let command = env::args()
        .nth(1)
        .unwrap_or_else(|| "dist-versioned".to_string());

    match command.as_str() {
        "dist-versioned" => {
            let dst = dist_versioned()?;
            println!("created {}", dst.display());
            Ok(())
        }
        "install-versioned" => install_versioned(),
        _ => {
            eprintln!(
                "usage: cargo run --release --bin xtask -- <dist-versioned|install-versioned>"
            );
            exit(2);
        }
    }
}

fn dist_versioned() -> Result<PathBuf, Box<dyn Error>> {
    build_yi_corr()?;

    let version = package_version()?;
    let src = PathBuf::from("target/release/yi-corr");
    let dst = PathBuf::from(format!("target/release/yi-corr-v{version}"));

    copy_executable(&src, &dst)?;
    Ok(dst)
}

fn install_versioned() -> Result<(), Box<dyn Error>> {
    let src = dist_versioned()?;
    let bin_dir = cargo_bin_dir()?;
    fs::create_dir_all(&bin_dir)?;

    let dst = bin_dir.join(src.file_name().ok_or("missing versioned binary name")?);
    copy_executable(&src, &dst)?;

    println!("installed {}", dst.display());
    Ok(())
}

fn build_yi_corr() -> Result<(), Box<dyn Error>> {
    let cargo = env::var_os("CARGO").unwrap_or_else(|| "cargo".into());
    let status = Command::new(cargo)
        .args(["build", "--release", "--offline", "--bin", "yi-corr"])
        .status()?;

    if !status.success() {
        return Err(format!("cargo build failed with {status}").into());
    }

    Ok(())
}

fn package_version() -> Result<String, Box<dyn Error>> {
    let manifest = fs::read_to_string("Cargo.toml")?;
    let mut in_package = false;
    for line in manifest.lines() {
        let trimmed = line.trim();
        if trimmed == "[package]" {
            in_package = true;
            continue;
        }
        if in_package && trimmed.starts_with('[') {
            break;
        }
        if in_package && trimmed.starts_with("version") {
            let (_, value) = trimmed
                .split_once('=')
                .ok_or("invalid package version line in Cargo.toml")?;
            return Ok(value.trim().trim_matches('"').to_string());
        }
    }
    Err("package version not found in Cargo.toml".into())
}

fn cargo_bin_dir() -> Result<PathBuf, Box<dyn Error>> {
    if let Some(cargo_home) = env::var_os("CARGO_HOME") {
        return Ok(PathBuf::from(cargo_home).join("bin"));
    }

    let home = env::var_os("HOME").ok_or("HOME is not set and CARGO_HOME is not set")?;
    Ok(PathBuf::from(home).join(".cargo/bin"))
}

fn copy_executable(src: &Path, dst: &Path) -> Result<(), Box<dyn Error>> {
    fs::copy(src, dst)?;

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;

        let mode = fs::metadata(src)?.permissions().mode();
        fs::set_permissions(dst, fs::Permissions::from_mode(mode))?;
    }

    Ok(())
}

use std::fs;
use std::path::Path;

fn parse_cache_size_bytes(text: &str) -> Option<u64> {
    let s = text.trim();
    if s.is_empty() {
        return None;
    }
    let last = s.chars().last()?;
    let (num, mul) = match last.to_ascii_uppercase() {
        'K' => (&s[..s.len() - 1], 1024u64),
        'M' => (&s[..s.len() - 1], 1024u64 * 1024u64),
        'G' => (&s[..s.len() - 1], 1024u64 * 1024u64 * 1024u64),
        _ => (s, 1u64),
    };
    num.trim().parse::<u64>().ok().map(|v| v.saturating_mul(mul))
}

fn detect_linux_l3_cache_bytes() -> Option<u64> {
    let base = Path::new("/sys/devices/system/cpu/cpu0/cache");
    let entries = fs::read_dir(base).ok()?;
    let mut best: Option<u64> = None;

    for entry in entries.flatten() {
        let p = entry.path();
        let level = match fs::read_to_string(p.join("level"))
            .ok()
            .and_then(|v| v.trim().parse::<u32>().ok())
        {
            Some(v) => v,
            None => continue,
        };
        if level < 3 {
            continue;
        }

        let cache_type = match fs::read_to_string(p.join("type")) {
            Ok(v) => v.trim().to_ascii_lowercase(),
            Err(_) => continue,
        };
        if cache_type != "unified" && cache_type != "data" {
            continue;
        }

        let size_bytes = match fs::read_to_string(p.join("size"))
            .ok()
            .and_then(|v| parse_cache_size_bytes(&v))
        {
            Some(v) if v > 0 => v,
            _ => continue,
        };

        best = Some(best.map_or(size_bytes, |cur| cur.max(size_bytes)));
    }

    best
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    let logical_cpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(0);
    let l3_bytes = detect_linux_l3_cache_bytes().unwrap_or(0);

    println!("cargo:rustc-env=YI_BUILD_LOGICAL_CPUS={logical_cpus}");
    println!("cargo:rustc-env=YI_BUILD_L3_CACHE_BYTES={l3_bytes}");
}

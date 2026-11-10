//! Detailed correlation diagnostics go to a ZIP; explicit status stays visible.
use std::fmt::Arguments;
use std::fs::{self, File};
use std::io::{self, BufWriter, IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use zip::{write::SimpleFileOptions, CompressionMethod, ZipWriter};

#[derive(Default)]
struct RuntimeLog {
    info: Option<ZipWriter<BufWriter<File>>>,
    console: Option<BufWriter<File>>,
    write_error: Option<io::Error>,
    terminal: bool,
    progress_line: bool,
    last_progress: Option<Instant>,
}

fn logger() -> &'static Mutex<RuntimeLog> {
    static LOGGER: OnceLock<Mutex<RuntimeLog>> = OnceLock::new();
    LOGGER.get_or_init(|| Mutex::new(RuntimeLog::default()))
}

impl RuntimeLog {
    fn record(&mut self, text: &str, newline: bool) {
        if let Some(writer) = self.info.as_mut() {
            let result = write_text(writer, text, newline);
            self.remember_error(result);
        }
    }

    fn remember_error(&mut self, result: io::Result<()>) {
        if let Err(error) = result {
            if self.write_error.is_none() {
                std::eprintln!("[warn] correlation log write failed: {error}");
                self.write_error = Some(error);
            }
        }
    }

    fn end_progress(&mut self) {
        if self.progress_line {
            let mut out = io::stdout().lock();
            let _ = out.write_all(b"\n");
            let _ = out.flush();
            self.progress_line = false;
        }
    }

    fn echo(&mut self, text: &str, newline: bool) {
        self.end_progress();
        let mut out = io::stdout().lock();
        let _ = write_text(&mut out, text, newline);
        let _ = out.flush();
        self.save_console(text, newline);
    }

    fn save_console(&mut self, text: &str, newline: bool) {
        if let Some(writer) = self.console.as_mut() {
            let result = write_text(writer, text, newline).and_then(|_| writer.flush());
            self.remember_error(result);
        }
    }
}

fn write_text(writer: &mut impl Write, text: &str, newline: bool) -> io::Result<()> {
    writer.write_all(text.as_bytes())?;
    if newline {
        writer.write_all(b"\n")?;
    }
    Ok(())
}

pub fn init(
    directory: &Path,
    save_stdout: bool,
) -> Result<(PathBuf, Option<PathBuf>), crate::utils::DynError> {
    fs::create_dir_all(directory)?;
    let now = SystemTime::now().duration_since(UNIX_EPOCH)?;
    let tag = crate::unix_seconds_to_yyyydddhhmmss(now.as_secs() as i64)?;
    let suffix = format!("{}_{:09}_{}", tag, now.subsec_nanos(), std::process::id());
    let info_path = directory.join(format!("yi-corr_{suffix}.info"));
    let mut info = ZipWriter::new(BufWriter::new(File::create_new(&info_path)?));
    // A single ordinary ZIP32 Deflate member works with both zipcat and zedit.
    info.start_file(
        "info.txt",
        SimpleFileOptions::default()
            .compression_method(CompressionMethod::Deflated)
            .compression_level(Some(1))
            .large_file(false),
    )?;
    writeln!(
        info,
        "# yi-corr {} details; ZIP/Deflate; started {} UTC",
        env!("CARGO_PKG_VERSION"),
        tag
    )?;
    writeln!(
        info,
        "# argv: {:?}",
        std::env::args_os().collect::<Vec<_>>()
    )?;
    writeln!(info, "# cwd: {}", std::env::current_dir()?.display())?;
    info.flush()?;
    let stdout_path =
        save_stdout.then(|| PathBuf::from("stdout").join(format!("stdout_{suffix}.log")));
    let console = if let Some(path) = &stdout_path {
        fs::create_dir_all(path.parent().unwrap())?;
        Some(BufWriter::new(File::create_new(path)?))
    } else {
        None
    };
    let mut log = logger()
        .lock()
        .map_err(|_| io::Error::other("correlation log lock poisoned"))?;
    *log = RuntimeLog {
        info: Some(info),
        console,
        terminal: io::stdout().is_terminal(),
        ..RuntimeLog::default()
    };
    Ok((info_path, stdout_path))
}

pub fn is_active() -> bool {
    logger()
        .lock()
        .map(|log| log.info.is_some())
        .unwrap_or(false)
}

pub fn detail(args: Arguments<'_>, newline: bool) {
    let text = args.to_string();
    if let Ok(mut log) = logger().lock() {
        log.record(&text, newline);
        let message = text.trim_start();
        if log.info.is_none() || message.starts_with("[warn]") || message.starts_with("[error]") {
            log.echo(&text, newline);
        }
    }
}

pub fn status(args: Arguments<'_>, newline: bool) {
    let text = args.to_string();
    if let Ok(mut log) = logger().lock() {
        log.record(&text, newline);
        log.echo(&text, newline);
    }
}

pub fn stderr(args: Arguments<'_>) {
    let text = args.to_string();
    if let Ok(mut log) = logger().lock() {
        log.record(&text, true);
        log.end_progress();
        let _ = write_text(&mut io::stderr().lock(), &text, true);
    }
}

pub fn progress(args: Arguments<'_>, complete: bool) {
    if let Ok(mut log) = logger().lock() {
        let interval = Duration::from_millis(if log.terminal { 500 } else { 10_000 });
        if !complete
            && log
                .last_progress
                .is_some_and(|last| last.elapsed() < interval)
        {
            return;
        }
        log.last_progress = (!complete).then(Instant::now);
        let text = args.to_string();
        log.record(&text, true);
        if log.terminal {
            let mut out = io::stdout().lock();
            let _ = write!(out, "\r{text}");
            if complete {
                let _ = out.write_all(b"\n");
            }
            let _ = out.flush();
            log.progress_line = !complete;
            log.save_console(&text, true);
        } else {
            log.echo(&text, true);
        }
        if let Some(info) = log.info.as_mut() {
            let result = info.flush();
            log.remember_error(result);
        }
    }
}

/// Always called after run_cli, including errors, to finalize the ZIP directory.
pub fn finish(error: Option<&str>) -> io::Result<()> {
    let mut log = logger()
        .lock()
        .map_err(|_| io::Error::other("correlation log lock poisoned"))?;
    log.end_progress();
    if let Some(error) = error {
        log.record(&format!("[error] {error}"), true);
    }
    log.record(
        if error.is_some() {
            "# status=error"
        } else {
            "# status=ok"
        },
        true,
    );
    if let Some(writer) = log.info.take() {
        let result = writer
            .finish()
            .map_err(io::Error::other)
            .and_then(|mut file| file.flush());
        log.remember_error(result);
    }
    if let Some(mut writer) = log.console.take() {
        let result = writer.flush();
        log.remember_error(result);
    }
    log.write_error.take().map_or(Ok(()), Err)
}

use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{mpsc, Arc};
use std::thread;
use std::time::Instant;

use crate::utils::DynError;
use crate::{compute_frame_delay_entry, DelayEvalConfig, FrameDelayEntry, PackedSampleReader};

pub(crate) struct Sector {
    pub frames: usize,
    pub starts: [u64; 2],
    pub samples: [u64; 2],
    pub d_seek: f64,
}

pub(crate) struct Block {
    pub sector: usize,
    pub frame: usize,
    pub starts: [u64; 2],
    // Original sector-relative frame origin minus this read window's origin.
    pub offsets: [i64; 2],
    pub raw: [Vec<u8>; 2],
    pub delays: Vec<FrameDelayEntry>,
    pub delay_s: f64,
}

impl Block {
    fn empty() -> Self {
        Self {
            sector: 0,
            frame: 0,
            starts: [0; 2],
            offsets: [0; 2],
            raw: [Vec::new(), Vec::new()],
            delays: Vec::new(),
            delay_s: 0.0,
        }
    }
}

// Use the exact per-frame integer shifts, not an assumed maximum delay drift.
// Clipping to the original sector preserves its existing boundary padding.
fn read_window(
    delays: &[FrameDelayEntry],
    frame: usize,
    fft: usize,
    samples: u64,
    word: u64,
    antenna: usize,
) -> (u64, u64, i64) {
    let base = frame as i128 * fft as i128;
    let mut lo = samples as i128;
    let mut hi = 0_i128;
    for (i, d) in delays.iter().enumerate() {
        let shift = if antenna == 0 { d.int1 } else { d.int2 };
        let begin = base + i as i128 * fft as i128 - shift as i128;
        let end = begin + fft as i128;
        if end > 0 && begin < samples as i128 {
            lo = lo.min(begin.max(0));
            hi = hi.max(end.min(samples as i128));
        }
    }
    if hi <= lo {
        return (0, 0, i64::try_from(base).expect("frame origin exceeds i64"));
    }
    let start = lo as u64 / word * word;
    let end = ((hi as u64).div_ceil(word) * word).min(samples);
    (
        start,
        end - start,
        i64::try_from(base - start as i128).expect("read offset exceeds i64"),
    )
}

fn read_into(
    reader: &mut PackedSampleReader,
    buf: &mut Vec<u8>,
    start: u64,
    samples: u64,
    bits: usize,
) -> Result<(), DynError> {
    let bytes = usize::try_from((samples as u128 * bits as u128).div_ceil(8))?;
    if bytes > buf.capacity() {
        buf.try_reserve_exact(bytes - buf.len())?;
    }
    buf.resize(bytes, 0);
    let start_bits = start
        .checked_mul(bits as u64)
        .ok_or("input bit offset overflow")?;
    reader.seek_to(start_bits / 8, (start_bits % 8) as u8)?;
    reader.read_packed_with_padding(buf)
}

pub(crate) struct Pipeline {
    ready: Option<mpsc::Receiver<Result<Block, String>>>,
    recycle: Option<mpsc::SyncSender<Block>>,
    handle: Option<thread::JoinHandle<()>>,
}

impl Pipeline {
    #[allow(clippy::too_many_arguments)]
    pub fn start(
        paths: [PathBuf; 2],
        bits: [usize; 2],
        fft: usize,
        sectors: Vec<Sector>,
        cfg: DelayEvalConfig,
        chunk_frames: usize,
        depth: usize,
        concurrent: bool,
        core: Option<core_affinity::CoreId>,
        produced: Arc<AtomicUsize>,
        produced_bytes: Arc<AtomicU64>,
    ) -> Self {
        let (ready_tx, ready) = mpsc::sync_channel(depth);
        // Fixed number of owned slots: queued, in computation, in the reader.
        let (recycle, free) = mpsc::sync_channel(depth + 2);
        for _ in 0..depth + 2 {
            recycle.send(Block::empty()).unwrap();
        }
        let handle = thread::spawn(move || {
            if let Some(core) = core {
                let _ = crate::affinity::set_current_thread_core(core);
            }
            let result = (|| -> Result<(), DynError> {
                let mut r1 = PackedSampleReader::open(&paths[0], 0, 0)?;
                let r2 = PackedSampleReader::open(&paths[1], 0, 0)?;
                thread::scope(|scope| -> Result<(), DynError> {
                    // Optional second persistent reader for separate USB devices.
                    let (request_tx, request_rx) = mpsc::sync_channel::<(Vec<u8>, u64, u64)>(0);
                    let (reply_tx, reply_rx) = mpsc::sync_channel::<Result<Vec<u8>, DynError>>(0);
                    let mut serial_r2 = Some(r2);
                    if concurrent {
                        let mut r2 = serial_r2.take().unwrap();
                        scope.spawn(move || {
                            while let Ok((mut buf, start, samples)) = request_rx.recv() {
                                let result = read_into(&mut r2, &mut buf, start, samples, bits[1])
                                    .map(|()| buf);
                                if reply_tx.send(result).is_err() {
                                    break;
                                }
                            }
                        });
                    }
                    let mut emitted = 0;
                    for (si, sector) in sectors.iter().enumerate() {
                        for frame in (0..sector.frames).step_by(chunk_frames) {
                            let Ok(mut block) = free.recv() else {
                                return Ok(());
                            };
                            let nf = chunk_frames.min(sector.frames - frame);
                            block.sector = si;
                            block.frame = frame;
                            let t = Instant::now();
                            block.delays.clear();
                            block.delays.extend((0..nf).map(|i| {
                                compute_frame_delay_entry(emitted + frame + i, &cfg, sector.d_seek)
                            }));
                            block.delay_s = t.elapsed().as_secs_f64();
                            let mut counts = [0; 2];
                            for ant in 0..2 {
                                let (start, count, offset) = read_window(
                                    &block.delays,
                                    frame,
                                    fft,
                                    sector.samples[ant],
                                    32 / bits[ant] as u64,
                                    ant,
                                );
                                block.starts[ant] = sector.starts[ant] + start;
                                block.offsets[ant] = offset;
                                counts[ant] = count;
                            }
                            if concurrent {
                                request_tx.send((
                                    std::mem::take(&mut block.raw[1]),
                                    block.starts[1],
                                    counts[1],
                                ))?;
                            }
                            // Always collect the concurrent result before returning an error,
                            // so the worker cannot remain blocked on its rendezvous send.
                            let first = read_into(
                                &mut r1,
                                &mut block.raw[0],
                                block.starts[0],
                                counts[0],
                                bits[0],
                            );
                            let second = if concurrent {
                                reply_rx.recv()?.map(|buf| {
                                    block.raw[1] = buf;
                                })
                            } else {
                                read_into(
                                    serial_r2.as_mut().unwrap(),
                                    &mut block.raw[1],
                                    block.starts[1],
                                    counts[1],
                                    bits[1],
                                )
                            };
                            first?;
                            second?;
                            let bytes = block.raw.iter().map(|v| v.len() as u64).sum();
                            if ready_tx.send(Ok(block)).is_err() {
                                return Ok(());
                            }
                            produced.fetch_add(1, Ordering::Relaxed);
                            produced_bytes.fetch_add(bytes, Ordering::Relaxed);
                        }
                        emitted += sector.frames;
                    }
                    Ok(())
                })
            })();
            if let Err(e) = result {
                let _ = ready_tx.send(Err(e.to_string()));
            }
        });
        Self {
            ready: Some(ready),
            recycle: Some(recycle),
            handle: Some(handle),
        }
    }

    pub fn recv(&self) -> Result<Block, DynError> {
        self.ready.as_ref().unwrap().recv()?.map_err(|e| e.into())
    }

    pub fn recycle(&self, block: Block) {
        // Once the final block has been read the producer can already have exited.
        let _ = self.recycle.as_ref().unwrap().send(block);
    }

    pub fn finish(mut self) -> Result<(), DynError> {
        self.ready.take();
        self.recycle.take();
        if self.handle.take().unwrap().join().is_err() {
            return Err("input pipeline reader panicked".into());
        }
        Ok(())
    }
}

impl Drop for Pipeline {
    fn drop(&mut self) {
        // Release both possible waits before joining, including early consumer errors.
        self.ready.take();
        self.recycle.take();
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::{build_decode_plan, decode_block_into_with_plan};
    use crate::{decode_shifted_frame_from_chunk, DecodeWindowScratch};

    #[test]
    fn chunk_windows_match_whole_sector_with_integer_delays_and_padding() {
        for bits in [1, 2, 4, 8] {
            let word = 32 / bits as u64;
            let fft = 32;
            let samples = 11 * fft + word as usize;
            let raw: Vec<u8> = (0..samples * bits / 8)
                .map(|i| (i * 53 + 79) as u8)
                .collect();
            let levels: Vec<f64> = (0..1usize << bits).map(|i| i as f64 - 1.5).collect();
            for map in [
                (0..32).collect::<Vec<_>>(),
                (0..32).map(|i| i ^ 7).collect(),
            ] {
                let plan = build_decode_plan(bits, &map, &levels).unwrap();
                for lsb in [false, true] {
                    for abs_start in [0, 1] {
                        let mut decoded = vec![0.; samples];
                        decode_block_into_with_plan(
                            &raw,
                            samples,
                            &plan,
                            &mut decoded,
                            lsb,
                            abs_start == 1,
                        )
                        .unwrap();
                        for frame in [0, 3, 6, 9] {
                            let nf = 3.min(11 - frame);
                            for shift in [-1000, -65, -17, -1, 0, 1, 17, 65, 1000] {
                                let delays: Vec<_> = (0..nf)
                                    .map(|i| FrameDelayEntry {
                                        int1: shift + i as i64 * 3,
                                        ..Default::default()
                                    })
                                    .collect();
                                let (start, count, offset) =
                                    read_window(&delays, frame, fft, samples as u64, word, 0);
                                let bytes = &raw[start as usize * bits / 8
                                    ..(start + count) as usize * bits / 8];
                                let mut scratch = DecodeWindowScratch::new();
                                let mut actual = vec![123.; fft];
                                for (i, d) in delays.iter().enumerate() {
                                    decode_shifted_frame_from_chunk(
                                        bytes,
                                        abs_start + start,
                                        i,
                                        fft,
                                        bits,
                                        word,
                                        &plan,
                                        lsb,
                                        d.int1 - offset,
                                        &mut actual,
                                        &mut scratch,
                                    )
                                    .unwrap();
                                    for k in 0..fft {
                                        let index = ((frame + i) * fft + k) as i64 - d.int1;
                                        let expected = if index >= 0 && index < samples as i64 {
                                            decoded[index as usize]
                                        } else {
                                            0.
                                        };
                                        assert_eq!(
                                            actual[k].to_bits(),
                                            expected.to_bits(),
                                            "bits={bits} frame={frame} i={i} shift={shift} k={k}"
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn cancellation_releases_blocked_readers_and_open_errors_propagate() {
        use std::time::{Duration, SystemTime, UNIX_EPOCH};
        let root = std::env::temp_dir().join(format!(
            "fx-input-drop-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&root).unwrap();
        let file = root.join("samples.raw");
        std::fs::write(&file, vec![0_u8; 4096]).unwrap();
        for concurrent in [false, true] {
            for consume in [false, true] {
                let paths = [file.clone(), file.clone()];
                let (done_tx, done_rx) = mpsc::channel();
                thread::spawn(move || {
                    let p = Pipeline::start(
                        paths,
                        [2, 2],
                        32,
                        vec![Sector {
                            frames: 100,
                            starts: [0; 2],
                            samples: [3216; 2],
                            d_seek: 0.,
                        }],
                        DelayEvalConfig {
                            fs: 8192.,
                            frame_dt: 32. / 8192.,
                            fx_integer_delay: true,
                            ..Default::default()
                        },
                        3,
                        1,
                        concurrent,
                        None,
                        Arc::new(AtomicUsize::new(0)),
                        Arc::new(AtomicU64::new(0)),
                    );
                    if consume {
                        let b = p.recv().unwrap();
                        p.recycle(b);
                    }
                    // Give the producer time to block on a full ready queue.
                    thread::sleep(Duration::from_millis(20));
                    drop(p);
                    done_tx.send(()).unwrap();
                });
                done_rx
                    .recv_timeout(Duration::from_secs(5))
                    .expect("reader cancellation deadlocked");
            }
            let p = Pipeline::start(
                [file.clone(), root.join("missing.raw")],
                [2, 2],
                32,
                vec![Sector {
                    frames: 1,
                    starts: [0; 2],
                    samples: [48; 2],
                    d_seek: 0.,
                }],
                DelayEvalConfig {
                    fs: 8192.,
                    ..Default::default()
                },
                3,
                1,
                concurrent,
                None,
                Arc::new(AtomicUsize::new(0)),
                Arc::new(AtomicU64::new(0)),
            );
            assert!(p.recv().is_err());
            p.finish().unwrap();
        }
        std::fs::remove_dir_all(root).unwrap();
    }
}

//! Limit simultaneous FFT jobs using completed input and compute windows.
//! This is a throughput headroom policy, not a measurement of disk hardware speed.
pub struct Balance {
    max: usize,
    active: usize,
    chunks: usize,
    fast_chunks: usize,
    supply: f64,
    compute: f64,
    wait: f64,
}

impl Balance {
    pub fn new(max: usize) -> Self {
        let max = max.max(1);
        Self {
            max,
            active: max,
            chunks: 0,
            fast_chunks: 0,
            supply: 0.0,
            compute: 0.0,
            wait: 0.0,
        }
    }

    pub fn active(&self) -> usize {
        self.active
    }

    // Input time includes delay preparation but excludes producer queue waits.
    // Keep ~30% compute headroom and sample repeatedly: early cached input is
    // not assumed to represent the remainder of the observation.
    pub fn observe(&mut self, supply: f64, compute: f64, wait: f64) -> Option<(usize, usize)> {
        if [supply, compute, wait]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return None;
        }
        // Recover within two completed chunks when computation is clearly
        // behind the new supply rate; do not wait for the slower downshift window.
        if self.active < self.max && compute > supply * 1.25 && wait < compute * 0.1 {
            self.fast_chunks += 1;
        } else {
            self.fast_chunks = 0;
        }
        if self.fast_chunks >= 2 {
            let old = self.active;
            self.active = old.saturating_mul(2).min(self.max);
            self.fast_chunks = 0;
            self.chunks = 0;
            self.supply = 0.0;
            self.compute = 0.0;
            self.wait = 0.0;
            return Some((old, self.active));
        }
        self.chunks += 1;
        self.supply += supply;
        self.compute += compute;
        self.wait += wait;
        if self.chunks < 8 || self.compute + self.wait < 0.5 {
            return None;
        }
        let ratio = self.compute / self.supply.max(1e-9);
        let wait_fraction = self.wait / (self.compute + self.wait).max(1e-9);
        let target =
            ((self.active as f64 * ratio / 0.7).ceil() as usize).clamp(2.min(self.max), self.max);
        let old = self.active;
        if ratio > 0.85 && target > old {
            // Respond quickly when a previously slow reader speeds up.
            self.active = target.min(old.saturating_mul(2));
        } else if wait_fraction > 0.25 && target < old {
            // Do not reduce merely because the producer once filled its queue.
            self.active = target.max(old.div_ceil(2));
        }
        self.chunks = 0;
        self.supply = 0.0;
        self.compute = 0.0;
        self.wait = 0.0;
        (old != self.active).then_some((old, self.active))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn follows_slow_input_then_recovers_when_input_speeds_up() {
        let mut b = Balance::new(14);
        // Work cost is 0.26 CPU-seconds per chunk; input takes 0.1 s.
        for _ in 0..40 {
            let c = 0.26 / b.active() as f64;
            b.observe(0.1, c, (0.1 - c).max(0.0));
        }
        assert_eq!(b.active(), 4);
        for _ in 0..4 {
            let c = 0.26 / b.active() as f64;
            b.observe(0.005, c, 0.0);
        }
        assert_eq!(b.active(), 14);
    }

    #[test]
    fn transient_and_queued_input_do_not_reduce_parallelism() {
        let mut b = Balance::new(14);
        b.observe(1.0, 0.02, 1.0);
        assert_eq!(b.active(), 14);
        let mut b = Balance::new(14);
        for _ in 0..100 {
            b.observe(0.1, 0.02, 0.0);
        }
        assert_eq!(b.active(), 14);
        b.observe(f64::NAN, 1.0, 1.0);
        assert_eq!(b.active(), 14);
    }

    #[test]
    fn respects_small_cpu_budgets() {
        for max in [1, 2] {
            let mut b = Balance::new(max);
            for _ in 0..100 {
                b.observe(0.1, 0.001, 0.1);
            }
            assert_eq!(b.active(), max);
        }
    }
}

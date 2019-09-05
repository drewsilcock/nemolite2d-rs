use std::time::{Duration, Instant};

pub struct Timer {
    pub region: &'static str,
    times: Vec<Duration>,
    current_time: Option<Instant>,
}

impl Timer {
    pub fn new(region: &'static str, expected_num_iterations: usize) -> Self {
        Timer {
            region,
            times: Vec::with_capacity(expected_num_iterations),
            current_time: None,
        }
    }
    pub fn start(&mut self) {
        assert!(
            self.current_time.is_none(),
            "You need to call stop() before calling start() again in region '{}'.",
            self.region
        );
        self.current_time = Some(Instant::now());
    }

    pub fn stop(&mut self) {
        let elapsed = self
            .current_time
            .unwrap_or_else(|| {
                panic!(
                    "You need to call start() before calling stop() in region '{}'.",
                    self.region
                )
            })
            .elapsed();
        self.times.push(elapsed);
        self.current_time = None;
    }

    pub fn count(&self) -> u32 {
        self.times.len() as u32
    }

    pub fn total(&self) -> Duration {
        self.times.iter().sum()
    }

    pub fn mean(&self) -> Option<Duration> {
        let count = self.count();

        match count {
            iters if iters > 0 => Some(self.total() / count as u32),
            _ => None,
        }
    }

    pub fn std_dev(&self) -> Option<Duration> {
        match (self.mean(), self.count()) {
            (Some(mean_value), num_iters) if num_iters > 0 => {
                let variance = self
                    .times
                    .iter()
                    .map(|value| {
                        let diff = mean_value.as_nanos() as i128 - (*value).as_nanos() as i128;

                        (diff * diff) as u128
                    })
                    .sum::<u128>() as f64
                    / f64::from(num_iters);

                Some(Duration::from_nanos(variance.sqrt() as u64))
            }
            _ => None,
        }
    }
}

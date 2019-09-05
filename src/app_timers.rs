use crate::timer::Timer;

pub struct AppTimers {
    pub step: Timer,
    pub continuity: Timer,
    pub momentum: Timer,
    pub boundary_conditions: Timer,
    pub next: Timer,
}

impl AppTimers {
    pub fn new(total_num_steps: usize) -> Self {
        AppTimers {
            step: Timer::new("Total Step", total_num_steps),
            continuity: Timer::new("Continuity", total_num_steps),
            momentum: Timer::new("Momentum", total_num_steps),
            boundary_conditions: Timer::new("Boundary Conditions", total_num_steps),
            next: Timer::new("Next", total_num_steps),
        }
    }

    pub fn generate_report(&self) -> String {
        let report = format!(
            "\n\
             -------------------------------------------------------------------------------------\n\
             Region                        Count               Total        Average        Std Err\n\
             -------------------------------------------------------------------------------------\n\
             {:<30}{:<10}{:>15}{:>15}{:>15}\n\
             {:<30}{:<10}{:>15}{:>15}{:>15}\n\
             {:<30}{:<10}{:>15}{:>15}{:>15}\n\
             {:<30}{:<10}{:>15}{:>15}{:>15}\n\
             {:<30}{:<10}{:>15}{:>15}{:>15}\
             ",
            self.step.region,
            self.step.count(),
            format!("{:3.5?}", self.step.total()),
            format!("{:3.5?}", self.step.mean().unwrap()),
            format!("{:3.5?}", self.step.std_dev().unwrap()),
            self.continuity.region,
            self.continuity.count(),
            format!("{:3.5?}", self.continuity.total()),
            format!("{:3.5?}", self.continuity.mean().unwrap()),
            format!("{:3.5?}", self.continuity.std_dev().unwrap()),
            self.momentum.region,
            self.momentum.count(),
            format!("{:3.5?}", self.momentum.total()),
            format!("{:3.5?}", self.momentum.mean().unwrap()),
            format!("{:3.5?}", self.momentum.std_dev().unwrap()),
            self.boundary_conditions.region,
            self.boundary_conditions.count(),
            format!("{:3.5?}", self.boundary_conditions.total()),
            format!("{:3.5?}", self.boundary_conditions.mean().unwrap()),
            format!("{:3.5?}", self.boundary_conditions.std_dev().unwrap()),
            self.next.region,
            self.next.count(),
            format!("{:3.5?}", self.next.total()),
            format!("{:3.5?}", self.next.mean().unwrap()),
            format!("{:3.5?}", self.next.std_dev().unwrap()),
        );

        report
    }

    pub fn generate_timings_csv(&self) -> String {
        let rows = [
            [
                "Region".to_owned(),
                "Count".to_owned(),
                "Total".to_owned(),
                "Average".to_owned(),
                "Std Dev".to_owned(),
            ],
            [
                self.step.region.to_owned(),
                format!("{}", self.step.count()),
                format!("{:?}", self.step.total()),
                format!("{:?}", self.step.mean().unwrap()),
                format!("{:?}", self.step.std_dev().unwrap()),
            ],
            [
                self.continuity.region.to_owned(),
                format!("{}", self.continuity.count()),
                format!("{:?}", self.continuity.total()),
                format!("{:?}", self.continuity.mean().unwrap()),
                format!("{:?}", self.continuity.std_dev().unwrap()),
            ],
            [
                self.momentum.region.to_owned(),
                format!("{}", self.momentum.count()),
                format!("{:?}", self.momentum.total()),
                format!("{:?}", self.momentum.mean().unwrap()),
                format!("{:?}", self.momentum.std_dev().unwrap()),
            ],
            [
                self.boundary_conditions.region.to_owned(),
                format!("{}", self.boundary_conditions.count()),
                format!("{:?}", self.boundary_conditions.total()),
                format!("{:?}", self.boundary_conditions.mean().unwrap()),
                format!("{:?}", self.boundary_conditions.std_dev().unwrap()),
            ],
            [
                self.next.region.to_owned(),
                format!("{}", self.next.count()),
                format!("{:?}", self.next.total()),
                format!("{:?}", self.next.mean().unwrap()),
                format!("{:?}", self.next.std_dev().unwrap()),
            ],
        ];

        let field_separater = ",";
        let line_separator = "\n";

        rows.iter()
            .map(|row| row.join(field_separater))
            .collect::<Vec<String>>()
            .join(line_separator)
    }
}

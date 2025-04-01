use std::collections::HashMap;
const IGRFCOEFFS: &str = include_str!("../../coeffs/shc/igrf14coeffs.txt");
const INTERVAL: f64 = 5.;
struct CoeffDetails {
    nmax: i16,
    coeffs: Vec<f64>,
}
pub(crate) struct IGRFCoeffs {
    coeffs: HashMap<i16, CoeffDetails>,
}

pub fn igrf_data() -> IGRFCoeffs {
    let parts = IGRFCOEFFS
        .split('\n')
        .filter(|s| !s.is_empty())
        .filter(|s| !s.starts_with('#'))
        .map(|s| s.split_whitespace().collect::<Vec<_>>())
        .collect::<Vec<_>>();

    let mut coeffs = parts[0]
        .iter()
        .skip(3) // skip to year columns
        .zip(parts[1].iter().skip(3))
        .enumerate()
        .map(|(i, (_, b))| {
            //let hdr = format!("{} {}", a, b);
            let epoch = if b.contains('-') {
                2000.
                    + b.split('-')
                        .last()
                        .and_then(|f| f.parse::<f32>().ok())
                        .unwrap()
            } else {
                b.parse::<f32>().unwrap_or_default()
            };
            let nmax = if epoch < 2000. {
                10
            } else if epoch > 2020. {
                8
            } else {
                13
            };

            (
                epoch as i16,
                CoeffDetails {
                    nmax,
                    coeffs: parts
                        .iter()
                        .skip(2)
                        .map(|line| line[i + 3].parse::<f64>().unwrap())
                        .collect::<Vec<_>>(),
                },
            )
        })
        .collect::<HashMap<i16, CoeffDetails>>();

    //todo: NASTY
    let bses = coeffs.get(&2025).unwrap().coeffs.clone();
    let ases = coeffs.get_mut(&2030).unwrap();

    ases.coeffs = ases
        .coeffs
        .iter()
        .zip(bses.iter())
        .map(|(a, b)| b + a * INTERVAL)
        .collect();

    IGRFCoeffs { coeffs }
}

fn find_date_factor(start_epoch: i16, end_epoch: i16, date: f64) -> f64 {
    fn secs_in_year(year: i16) -> i32 {
        let is_leap = year % 400 == 0 || (year % 4 == 0 && year % 100 != 0);
        if is_leap {
            366 * 24 * 3600
        } else {
            365 * 24 * 3600
        }
    }

    let dte2 = end_epoch as f64;
    if end_epoch <= start_epoch {
        0.
    } else if date > dte2 {
        let dte1 = start_epoch as f64;
        (date - dte1) / (dte2 - dte1)
    } else {
        let loc_interval = end_epoch - start_epoch;
        let mut total_secs = 0f64;
        let mut fraction_secs = 0f64;
        for x in 0..loc_interval {
            let year = start_epoch + x;
            let secs_in_year = secs_in_year(year) as f64;
            if year == date as i16 {
                let fraction_coeff = date - date.floor();
                fraction_secs = total_secs + fraction_coeff * secs_in_year;
            }
            total_secs += secs_in_year;
        }
        fraction_secs / total_secs
    }
}

impl IGRFCoeffs {
    fn extrapolate_coeffs(&self, start_epoch: i16, end_epoch: i16, date: f64) -> Vec<f64> {
        let start = self.coeffs.get(&start_epoch).unwrap();
        let end = self.coeffs.get(&end_epoch).unwrap();
        if start.nmax <= end.nmax {
            panic!("nmax1 <= nmax2")
        }
        let k = end.nmax * (end.nmax + 2);
        let l = start.nmax * (start.nmax + 2);
        start
            .coeffs
            .iter()
            .zip(end.coeffs.iter())
            .enumerate()
            .map(|(i, (coeff_start, coeff_end))| {
                if (i as i16) >= k && (i as i16) < l {
                    *coeff_start
                } else {
                    let sv = (coeff_end - coeff_start) / INTERVAL;
                    let factor = date - (start_epoch as f64);
                    coeff_start + factor * sv
                }
            })
            .collect::<Vec<f64>>()
    }

    fn interpolate_coeffs(&self, start_epoch: i16, end_epoch: i16, date: f64) -> (Vec<f64>, i16) {
        let factor = find_date_factor(start_epoch, end_epoch, date);
        let start = self.coeffs.get(&start_epoch).unwrap();
        let end = self.coeffs.get(&end_epoch).unwrap();

        let (nmax, k, l, interp) = match start.nmax.cmp(&end.nmax) {
            //before 2000.0
            std::cmp::Ordering::Equal => (start.nmax, start.nmax * (start.nmax + 2), -100, true),

            // between 1995.0 and 2000.0
            std::cmp::Ordering::Less => (
                end.nmax,
                start.nmax * (start.nmax + 2),
                end.nmax * (end.nmax + 2),
                false,
            ),

            // the last column has degree of 8
            // now it's anything after 2020.0
            std::cmp::Ordering::Greater => (
                start.nmax,
                end.nmax * (end.nmax + 2),
                start.nmax * (start.nmax + 2),
                true,
            ),
        };

        let values = start
            .coeffs
            .iter()
            .zip(end.coeffs.iter())
            .enumerate()
            .map(|(i, (coeff_start, coeff_end))| {
                if k <= (i as i16) && (i as i16) < l {
                    if interp {
                        *coeff_start
                    } else {
                        factor * coeff_end
                    }
                } else {
                    coeff_start + factor * (coeff_end - coeff_start)
                }
            })
            .collect::<Vec<_>>();
        (values, nmax)
    }

    pub(crate) fn coeffs(&self, date: f64) -> (Vec<f64>, Vec<f64>, i16) {
        let years = self.coeffs.keys().collect::<Vec<_>>();
        if !(1900. ..=2030.).contains(&date) {
            panic!("Date out of range");
        }

        let end = **years
            .iter()
            .filter(|x| (***x as f64) >= date)
            .min()
            .unwrap();

        let start = **years
            .iter()
            .filter(|x| (***x as f64) <= date)
            .max()
            .unwrap();

        let mac_epoch = **years.iter().max().unwrap();

        if date > 2030. {
            todo!("2030+")
        }

        let (coeffs_start, nmax) = self.interpolate_coeffs(start, end, date);

        if date + 1. < mac_epoch as f64 {
            let (coeffs_end, nmax) = self.interpolate_coeffs(start, end, date + 1.);
            (coeffs_start, coeffs_end, nmax)
        } else {
            let coeffs_end = self.extrapolate_coeffs(start, end, date + 1.);
            (coeffs_start, coeffs_end, nmax)
        }
    }
}

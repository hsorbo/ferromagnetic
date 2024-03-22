use std::vec;
const EARTH_SEMI_MAJOR_AXIS_KM: f32 = 6378.137;
const EARTH_SEMI_MINOR_AXIS_KM: f32 = 6356.7524;
const EARTH_REFERENCE_RADIUS_KM: f32 = 6371.2;
const G_COEFF: [[f32; 13]; 13] = [
    [
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        -29438.5, -1501.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        -2445.3, 3012.5, 1676.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        1351.1, -2352.3, 1225.6, 581.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        907.2, 813.7, 120.3, -335.0, 70.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        -232.6, 360.1, 192.4, -141.0, -157.4, 4.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        69.5, 67.4, 72.8, -129.8, -29.0, 13.2, -70.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        81.6, -76.1, -6.8, 51.9, 15.0, 9.3, -2.8, 6.7, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        24.0, 8.6, -16.9, -3.2, -20.6, 13.3, 11.7, -16.0, -2.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        5.4, 8.8, 3.1, -3.1, 0.6, -13.3, -0.1, 8.7, -9.1, -10.5, 0.0, 0.0, 0.0,
    ],
    [
        -1.9, -6.5, 0.2, 0.6, -0.6, 1.7, -0.7, 2.1, 2.3, -1.8, -3.6, 0.0, 0.0,
    ],
    [
        3.1, -1.5, -2.3, 2.1, -0.9, 0.6, -0.7, 0.2, 1.7, -0.2, 0.4, 3.5, 0.0,
    ],
    [
        -2.0, -0.3, 0.4, 1.3, -0.9, 0.9, 0.1, 0.5, -0.4, -0.4, 0.2, -0.9, 0.0,
    ],
];
const H_COEFF: [[f32; 13]; 13] = [
    [
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 4796.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, -2845.6, -642.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, -115.3, 245.0, -538.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 283.4, -188.6, 180.9, -329.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 47.4, 196.9, -119.4, 16.1, 100.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, -20.7, 33.2, 58.8, -66.5, 7.3, 62.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, -54.1, -19.4, 5.6, 24.4, 3.3, -27.5, -2.3, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 10.2, -18.1, 13.2, -14.6, 16.2, 5.7, -9.1, 2.2, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, -21.6, 10.8, 11.7, -6.8, -6.9, 7.8, 1.0, -3.9, 8.5, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 3.3, -0.3, 4.6, 4.4, -7.9, -0.6, -4.1, -2.8, -1.1, -8.7, 0.0, 0.0,
    ],
    [
        0.0, -0.1, 2.1, -0.7, -1.1, 0.7, -0.2, -2.1, -1.5, -2.5, -2.0, -2.3, 0.0,
    ],
    [
        0.0, -1.0, 0.5, 1.8, -2.2, 0.3, 0.7, -0.1, 0.3, 0.2, -0.9, -0.2, 0.7,
    ],
];

const DELTA_G: [[f32; 13]; 13] = [
    [
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        10.7, 17.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        -8.6, -3.3, 2.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        3.1, -6.2, -0.4, -10.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        -0.4, 0.8, -9.2, 4.0, -4.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        -0.2, 0.1, -1.4, 0.0, 1.3, 3.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        -0.5, -0.2, -0.6, 2.4, -1.1, 0.3, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.2, -0.2, -0.4, 1.3, 0.2, -0.4, -0.9, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 0.1, -0.5, 0.5, -0.2, 0.4, 0.2, -0.4, 0.3, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, -0.1, -0.1, 0.4, -0.5, -0.2, 0.1, 0.0, -0.2, -0.1, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 0.0, -0.1, 0.3, -0.1, -0.1, -0.1, 0.0, -0.2, -0.1, -0.2, 0.0, 0.0,
    ],
    [
        0.0, 0.0, -0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1, -0.1, 0.0,
    ],
    [
        0.1, 0.0, 0.0, 0.1, -0.1, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
];

const DELTA_H: [[f32; 13]; 13] = [
    [
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, -26.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, -27.1, -13.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 8.4, -0.4, 2.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, -0.6, 5.3, 3.0, -5.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 0.4, 1.6, -1.1, 3.3, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 0.0, -2.2, -0.7, 0.1, 1.0, 1.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 0.7, 0.5, -0.2, -0.1, -0.7, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, -0.3, 0.3, 0.3, 0.6, -0.1, -0.2, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0, -0.2, -0.1, -0.2, 0.1, 0.1, 0.0, -0.2, 0.4, 0.3, 0.0, 0.0, 0.0,
    ],
    [
        0.0, 0.1, -0.1, 0.0, 0.0, -0.2, 0.1, -0.1, -0.2, 0.1, -0.1, 0.0, 0.0,
    ],
    [
        0.0, 0.0, 0.1, 0.0, 0.1, 0.0, 0.0, 0.1, 0.0, -0.1, 0.0, -0.1, 0.0,
    ],
    [
        0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
];

const BASE_TIME: i64 = 1422745200000; // 1422748800000;

pub fn compute_schmidt_quasi_norm_factors(max_n: usize) -> Vec<Vec<f32>> {
    let mut schmidt_quasi_norm: Vec<Vec<f32>> = vec![];
    schmidt_quasi_norm.push(vec![1.0]);
    for n in 1..=max_n {
        let mut current = vec![];
        current.push(schmidt_quasi_norm[n - 1][0] * (2 * n - 1) as f32 / n as f32);
        for m in 1..=n {
            let a = (((n - m + 1) as f32 * if m == 1 { 2.0 } else { 1.0 }) / (n + m) as f32).sqrt();
            current.push(current.last().unwrap() * a);
        }
        schmidt_quasi_norm.push(current);
    }
    schmidt_quasi_norm
}

//const SCHMIDT_QUASI_NORM_FACTORS: Vec<Vec<f32>> = compute_schmidt_quasi_norm_factors(G_COEFF.len());

fn compute_geocentric_coordinates(
    gd_lat_deg: f32,
    gd_lon_deg: f32,
    alt_meters: f32,
) -> (f32, f32, f32) {
    let altitude_km = alt_meters / 1000.0;
    //let a2 = 4.0680636E7;
    let a2 = EARTH_SEMI_MAJOR_AXIS_KM * EARTH_SEMI_MAJOR_AXIS_KM;
    let b2 = EARTH_SEMI_MINOR_AXIS_KM * EARTH_SEMI_MINOR_AXIS_KM;
    //let b2 = 4.04083E7;
    let gd_lat_rad = gd_lat_deg.to_radians();
    let clat = gd_lat_rad.cos();
    let slat = gd_lat_rad.sin();
    let tlat = slat / clat;
    let lat_rad = (a2 * clat * clat + b2 * slat * slat).sqrt();
    let gc_lat_rad = (tlat * (lat_rad * altitude_km + b2) / (lat_rad * altitude_km + a2)).atan();
    let gc_lon_rad = gd_lon_deg.to_radians();
    let rad_sq = altitude_km * altitude_km
        + 2.0 * altitude_km * (a2 * clat * clat + b2 * slat * slat).sqrt()
        + (a2 * a2 * clat * clat + b2 * b2 * slat * slat) / (a2 * clat * clat + b2 * slat * slat);
    let gc_radius_km = rad_sq.sqrt();
    (gc_lat_rad, gc_lon_rad, gc_radius_km)
}

struct LegendreTable {
    p: Vec<Vec<f32>>,
    p_deriv: Vec<Vec<f32>>,
}

impl LegendreTable {
    pub fn new(max_n: usize, theta_rad: f32) -> Self {
        let cos = theta_rad.cos();
        let sin = theta_rad.sin();
        let mut t = LegendreTable {
            p: vec![vec![]; max_n + 1],
            p_deriv: vec![vec![]; max_n + 1],
        };
        t.p[0].push(1.0);
        t.p_deriv[0].push(0.0);
        for n in 1..=max_n {
            let mut current = vec![];
            let mut current_deriv = vec![];
            for m in 0..=n {
                if n == m {
                    current.push(sin * t.p[n - 1][m - 1]);
                    current_deriv.push(cos * t.p[n - 1][m - 1] + sin * t.p_deriv[n - 1][m - 1]);
                } else if n == 1 || m == n - 1 {
                    current.push(cos * t.p[n - 1][m]);
                    current_deriv.push(-sin * t.p[n - 1][m] + cos * t.p_deriv[n - 1][m]);
                } else {
                    assert!(n > 1 && m < n - 1);
                    let nf = n as f32;
                    let mf = m as f32;
                    let k = ((nf - 1.) * (nf - 1.) - mf * mf) / ((2. * nf - 1.) * (2. * nf - 3.));
                    current.push(cos * t.p[n - 1][m] - k * t.p[n - 2][m]);
                    current_deriv.push(
                        -sin * t.p[n - 1][m] + cos * t.p_deriv[n - 1][m] - k * t.p_deriv[n - 2][m],
                    );
                }
            }
            t.p[n] = current;
            t.p_deriv[n] = current_deriv;
        }
        t
    }
}

pub struct GeomagneticField {
    x: f32,
    y: f32,
    z: f32,
    gc_lat_rad: f32,
    gc_lon_rad: f32,
    gc_radius_km: f32,
}

impl GeomagneticField {
    pub fn new(gd_lat_deg: f32, gd_lon_deg: f32, alt_meters: f32, time_ms: i64) -> Self {
        //todo static
        let SCHMIDT_QUASI_NORM_FACTORS = compute_schmidt_quasi_norm_factors(G_COEFF.len());
        let max_n = G_COEFF.len();

        let mut field = {
            let lat_cut = 89.99999f32.min(gd_lat_deg).max(-89.99999);
            let (gc_lat_rad, gc_lon_rad, gc_radius_km) =
                compute_geocentric_coordinates(lat_cut, gd_lon_deg, alt_meters);
            Self {
                x: 0.0,
                y: 0.0,
                z: 0.0,
                gc_lat_rad: gc_lat_rad,
                gc_lon_rad: gc_lon_rad,
                gc_radius_km: gc_radius_km,
            }
        };

        assert!(G_COEFF.len() == H_COEFF.len());
        //let legendre = LegendreTable::new(max_n - 1, 1.5707963267948966 - kake.gc_lat_rad);
        let legendre =
            LegendreTable::new(max_n - 1, std::f32::consts::FRAC_PI_2 - field.gc_lat_rad);
        let mut relative_radius_power = vec![0.0; max_n + 2];
        relative_radius_power[0] = 1.0;
        relative_radius_power[1] = EARTH_REFERENCE_RADIUS_KM / field.gc_radius_km;
        for n in 2..=max_n + 1 {
            relative_radius_power[n] = relative_radius_power[n - 1] * relative_radius_power[1];
        }

        let mut sin_lon: Vec<f32> = vec![0.; max_n];
        let mut cos_lon: Vec<f32> = vec![0.; max_n];
        sin_lon[0] = 0.;
        cos_lon[0] = 1.;
        sin_lon[1] = field.gc_lon_rad.sin();
        cos_lon[1] = field.gc_lon_rad.cos();
        for m in 2..=max_n - 1 {
            let x = m >> 1;
            sin_lon[m] = sin_lon[m - x] * cos_lon[x] + cos_lon[m - x] * sin_lon[x];
            cos_lon[m] = cos_lon[m - x] * cos_lon[x] - sin_lon[m - x] * sin_lon[x];
        }
        let inverse_cos_latitude = 1.0 / field.gc_lat_rad.cos();
        let time_delta = time_ms - BASE_TIME;
        let years_since_base = (time_delta as f32) / 3.1536001E10; //(1000.0 * 60.0 * 60.0 * 24.0 * 365.25);
        let mut x = 0f32;
        let mut z = 0f32;
        for n in 1..max_n {
            for m in 0..=n {
                let g = G_COEFF[n][m] + DELTA_G[n][m] * years_since_base;
                let h = H_COEFF[n][m] + DELTA_H[n][m] * years_since_base;
                x += relative_radius_power[n + 2]
                    * (g * cos_lon[m] + h * sin_lon[m])
                    * legendre.p_deriv[n][m]
                    * SCHMIDT_QUASI_NORM_FACTORS[n][m];

                field.y += relative_radius_power[n + 2]
                    * (m as f32)
                    * (g * sin_lon[m] - h * cos_lon[m])
                    * legendre.p[n][m]
                    * SCHMIDT_QUASI_NORM_FACTORS[n][m]
                    * inverse_cos_latitude;

                z -= (n as f32 + 1.)
                    * relative_radius_power[n + 2]
                    * (g * cos_lon[m] + h * sin_lon[m])
                    * legendre.p[n][m]
                    * SCHMIDT_QUASI_NORM_FACTORS[n][m];
            }
            let lat_diff_rad = gd_lat_deg.to_radians() - field.gc_lat_rad;
            field.x = x * lat_diff_rad.cos() + z * lat_diff_rad.sin();
            field.z = -x * lat_diff_rad.sin() + z * lat_diff_rad.cos();
        }

        field
    }
    pub fn get_declination(&self) -> f32 {
        (self.y).atan2(self.x).to_degrees()
    }
}


#[cfg(test)]
mod tests {

    use crate::wmm::wmm_2015::compute_schmidt_quasi_norm_factors;

    use super::GeomagneticField;

    #[test]
    fn test() {
        let x = GeomagneticField::new(20.200628, -87.500114, 0.0, 1609545600000);
        let z = x.get_declination();
        assert_eq!(z, -1.7856733);
    }

    #[test]
    fn test2() {
        let x = compute_schmidt_quasi_norm_factors(6);
        assert_eq!(
            x,
            vec![
                vec![1.0],
                vec![1.0, 1.0],
                vec![1.5, 1.7320508, 0.8660254],
                vec![2.5, 3.0618622, 1.9364917, 0.7905695],
                vec![4.375, 5.533986, 3.913119, 2.09165, 0.73950994],
                vec![7.875, 10.166581, 7.6852126, 4.7062125, 2.21853, 0.70156074],
                vec![14.4375, 18.903124, 14.944232, 9.962822, 5.4568624, 2.326814, 0.6716933]
            ]
        );
    }
}

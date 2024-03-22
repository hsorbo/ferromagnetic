mod coeffs;
mod math;

#[derive(Debug)]
pub struct IGRFresults {
    pub declination: f64,
    pub declination_sv: f64,
    pub inclination: f64,
    pub inclination_sv: f64,
    pub horizontal_intensity: f64,
    pub horizontal_sv: f64,
    pub north_component: f64,
    pub north_sv: f64,
    pub east_component: f64,
    pub east_sv: f64,
    pub vertical_component: f64,
    pub vertical_sv: f64,
    pub total_intensity: f64,
    pub total_sv: f64,
}

pub struct IGRF {
    coeffs: coeffs::IGRF13Coeffs,
}
impl Default for IGRF {
    fn default() -> IGRF {
        IGRF {
            coeffs: coeffs::igrf_data(),
        }
    }
}
impl IGRF {
    pub fn calc(&self, lat: f64, lon: f64, alt: f64, date: f64) -> IGRFresults {
        //input validation
        let (start_coeffs, end_coeffs, nmax) = self.coeffs.coeffs(date);
        let (a, b) = math::shval3(lat, lon, alt, nmax as usize, &start_coeffs, &end_coeffs);
        let mut dif_a = math::Difh::default();
        let mut dif_b = math::Difh::default();
        dif_a.xyz(&a);
        dif_b.xyz(&b);

        let ddot = {
            let mut ddot = (dif_b.declination - dif_a.declination).to_degrees();
            if ddot > 180.0 {
                ddot -= 360.0
            }
            if ddot <= -180.0 {
                ddot += 360.0
            }
            ddot * 60.0
        };

        IGRFresults {
            declination: dif_a.declination.to_degrees(),
            declination_sv: ddot,
            inclination: dif_a.inclination.to_degrees(),
            inclination_sv: (dif_b.inclination - dif_a.inclination).to_degrees() * 60.,
            horizontal_intensity: dif_a.horizontal_intensity,
            horizontal_sv: dif_b.horizontal_intensity - dif_a.horizontal_intensity,
            north_component: a.x,
            north_sv: b.x - a.x,
            east_component: a.y,
            east_sv: b.y - a.y,
            vertical_component: a.z,
            vertical_sv: b.z - a.z,
            total_intensity: dif_a.total_intensity,
            total_sv: dif_b.total_intensity - dif_a.total_intensity,
        }
    }
}

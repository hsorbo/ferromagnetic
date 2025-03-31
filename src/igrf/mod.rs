use crate::{MagneticComponents, OrthogonalStrength};

mod coeffs;
mod math;

pub struct IGRFresults {
    pub result: MagneticComponents,
    // Annual changes
    pub sv: MagneticComponents,
}

pub struct IGRF {
    coeffs: coeffs::IGRFCoeffs,
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
        let dif_a = math::Difh::from_orthognal_strength(&a);
        let dif_b = math::Difh::from_orthognal_strength(&b);

        let result = MagneticComponents {
            declination: dif_a.declination.to_degrees(),
            inclination: dif_a.inclination.to_degrees(),
            horizontal_intensity: dif_a.horizontal_intensity,
            orthogonal_strength: a,
            total_intensity: dif_a.total_intensity,
        };

        let sv = MagneticComponents {
            declination: {
                let mut ddot = (dif_b.declination - dif_a.declination).to_degrees();
                if ddot > 180.0 {
                    ddot -= 360.0
                }
                if ddot <= -180.0 {
                    ddot += 360.0
                }
                ddot * 60.0
            },
            inclination: (dif_b.inclination - dif_a.inclination).to_degrees() * 60.,
            horizontal_intensity: dif_b.horizontal_intensity - dif_a.horizontal_intensity,
            orthogonal_strength: OrthogonalStrength {
                north: b.north - result.orthogonal_strength.north,
                east: b.east - result.orthogonal_strength.east,
                down: b.down - result.orthogonal_strength.down,
            },
            total_intensity: dif_b.total_intensity - dif_a.total_intensity,
        };

        IGRFresults { result, sv }
    }
}

use crate::OrthogonalStrength;

pub struct Difh {
    pub declination: f64,
    pub inclination: f64,
    pub horizontal_intensity: f64,
    pub total_intensity: f64,
}

impl Difh {
    pub fn from_orthognal_strength(p: &OrthogonalStrength) -> Self {
        const SN: f64 = 0.0001;

        let horizontal_intensity = (p.north * p.north + p.east * p.east).sqrt();
        let total_intensity = (p.north * p.north + p.east * p.east + p.down * p.down).sqrt();

        let inclination = if total_intensity < SN {
            std::f64::NAN
        } else {
            p.down.atan2(horizontal_intensity)
        };

        let declination = if total_intensity < SN || horizontal_intensity < SN {
            std::f64::NAN
        } else if horizontal_intensity + p.north < SN {
            std::f64::consts::PI
        } else {
            2.0 * p.east.atan2(horizontal_intensity + p.north)
        };
        Self {
            declination,
            inclination,
            horizontal_intensity,
            total_intensity,
        }
    }
}

//https://github.com/igp-gravity/geoist/blob/658aadab8074bffcbc6b3861671d35b3012502e9/geoist/pfm/igrf.py
//https://github.com/SuperDARN/rst/tree/734542037d86b24475d8da571ed668659a858e27/codebase/analysis/src.lib/igrf/igrf.1.13
//https://github.com/wallscavesurvey/walls/blob/master/geomag70/geomag70_org.c
//https://github.com/proway2/go-igrf

// // Computes field components from spherical harmonic (sh) models.
// // The calculation is performed for two sets of coeffs for a single location,
// // thus it returns two sets of X, Y, Z.
// //
// // X - northward component
// //
// // Y - eastward component
// //
// // Z - vertically-downward component
pub fn shval3(
    flat: f64,
    flon: f64,
    elev: f64,
    nmax: usize,
    gha: &[f64],
    ghb: &[f64],
) -> (OrthogonalStrength, OrthogonalStrength) {
    // 	// similar to shval3 from C implementation
    let earths_radius: f64 = 6371.2;
    let dtr: f64 = 0.01745329;
    // a2,b2     - squares of semi-major and semi-minor axes of
    // the reference spheroid used for transforming
    // between geodetic and geocentric coordinates or components
    let a2: f64 = 40680631.59; /* WGS84 */
    let b2: f64 = 40408299.98; /* WGS84 */

    let mut slat: f64 = (flat * dtr).sin();

    let mut clat = {
        let aa = if (90.0 - flat) < 0.001 {
            89.999
        } else if (90.0 + flat) < 0.001 {
            -89.999
        } else {
            flat
        };
        (aa * dtr).cos()
    };

    let mut sl: [f64; 14] = [0.0; 14];
    let mut cl: [f64; 14] = [0.0; 14];
    sl[1] = (flon * dtr).sin();
    cl[1] = (flon * dtr).cos();

    // this block is for geodetic coordinate system ->
    let aa = a2 * clat * clat;
    let bb = b2 * slat * slat;
    let cc = aa + bb;
    let dd = cc.sqrt();
    let r = (elev * (elev + 2.0 * dd) + (a2 * aa + b2 * bb) / cc).sqrt();
    let cd = (elev + dd) / r;
    let sd = (a2 - b2) / dd * slat * clat / r;
    {
        let old_slat = slat;
        slat = slat * cd - clat * sd;
        clat = clat * cd + old_slat * sd;
    }

    // <- this block is for geodetic coordinate system
    let ratio: f64 = earths_radius / r;

    let mut p: [f64; 119] = [0.0; 119];
    p[1] = 2.0 * slat;
    p[2] = 2.0 * clat;
    p[3] = 4.5 * slat * slat - 1.5;
    p[4] = 3.0 * 3.0f64.sqrt() * clat * slat;

    let mut q: [f64; 119] = [0.0; 119];
    q[1] = -clat;
    q[2] = slat;
    q[3] = -3.0 * clat * slat;
    q[4] = 3.0f64.sqrt() * (slat * slat - clat * clat);

    let mut l: usize = 0;
    let mut n: usize = 0;
    let mut m: usize = 1;
    let npq = (nmax * (nmax + 3)) / 2;

    let mut p2 = OrthogonalStrength::default();
    let mut p1 = OrthogonalStrength::default();

    let mut rr: f64 = 0.0;
    let mut fnn: f64 = 0.0;

    for k in 1..npq + 1 {
        if n < m {
            m = 0;
            n += 1;
            rr = ratio.powf((n + 2) as f64);
            fnn = n as f64;
        }
        let fm = m as f64;
        if k >= 5 {
            if m == n {
                let aa = (1.0 - 0.5 / fm).sqrt();
                let j = k - n - 1;
                p[k] = (1.0 + 1.0 / fm) * aa * clat * p[j];
                q[k] = aa * (clat * q[j] + slat / fm * p[j]);
                sl[m] = sl[m - 1] * cl[1] + cl[m - 1] * sl[1];
                cl[m] = cl[m - 1] * cl[1] - sl[m - 1] * sl[1];
            } else {
                let aa = (fnn * fnn - fm * fm).sqrt();
                let bb = ((fnn - 1.0) * (fnn - 1.0) - (fm * fm)).sqrt() / aa;
                let cc = (2.0 * fnn - 1.0) / aa;
                let ii = k - n;
                let j = k - 2 * n + 1;
                p[k] = (fnn + 1.0) * (cc * slat / fnn * p[ii] - bb / (fnn - 1.0) * p[j]);
                q[k] = cc * (slat * q[ii] - clat / fnn * p[ii]) - bb * q[j];
            }
        }

        let gh = &gha;
        let pp = &mut p1;
        if m == 0 {
            pp.north += (rr * gh[l]) * q[k];
            pp.down -= (rr * gh[l]) * p[k];
        } else {
            let b = rr * gh[l + 1];
            let c = (rr * gh[l]) * cl[m] + b * sl[m];
            pp.north += c * q[k];
            pp.down -= c * p[k];
            pp.east += if clat > 0.0 {
                ((rr * gh[l]) * sl[m] - b * cl[m]) * fm * p[k] / ((fnn + 1.0) * clat)
            } else {
                ((rr * gh[l]) * sl[m] - b * cl[m]) * q[k] * slat
            };
        }

        let gh = &ghb;
        let pp = &mut p2;
        if m == 0 {
            pp.north += (rr * gh[l]) * q[k];
            pp.down -= (rr * gh[l]) * p[k];
        } else {
            let b = rr * gh[l + 1];
            let c = (rr * gh[l]) * cl[m] + b * sl[m];
            pp.north += c * q[k];
            pp.down -= c * p[k];
            pp.east += if clat > 0.0 {
                ((rr * gh[l]) * sl[m] - b * cl[m]) * fm * p[k] / ((fnn + 1.0) * clat)
            } else {
                ((rr * gh[l]) * sl[m] - b * cl[m]) * q[k] * slat
            };
        }

        l += if m == 0 { 1 } else { 2 };
        m += 1;
    }
    {
        let old_x = p1.north;
        p1.north = p1.north * cd + p1.down * sd;
        p1.down = p1.down * cd - old_x * sd;
    }
    {
        let old_x = p2.north;
        p2.north = p2.north * cd + p2.down * sd;
        p2.down = p2.down * cd - old_x * sd;
    }

    (p1, p2)
}

#[cfg(test)]
mod tests {
    use float_eq::assert_float_eq;

    use super::*;

    #[test]
    fn difh_xyz() {
        let difh = Difh::from_orthognal_strength(&OrthogonalStrength {
            north: 7074.026894642207,
            east: 4596.188334222297,
            down: 62650.49087477549,
        });
        assert_float_eq!(difh.declination, 0.5761834859773884, rel <= 1e-12);
        assert_float_eq!(difh.inclination, 1.4369489460941476, rel <= 1e-12);
        assert_float_eq!(difh.horizontal_intensity, 8436.041945709041, rel <= 1e-12);
        assert_float_eq!(difh.total_intensity, 63215.90630972626, rel <= 1e-12);
    }

    #[test]
    fn test_shval3() {
        let start_coeffs: [f64; 195] = [
            -31535.104326396497,
            -2298.,
            5920.7007119386635,
            -682.0972070098576,
            2907.2987404162104,
            -1063.498630887185,
            935.6935925520263,
            1115.4030668127054,
            1023.499178532311,
            -1471.498630887185,
            -332.6985213581599,
            1254.3009309967142,
            6.09830230010953,
            578.2965498357065,
            518.7023548740416,
            876.3997809419496,
            629.499178532311,
            195.79956188389923,
            659.3003833515882,
            -69.79956188389923,
            -362.89895947426066,
            -209.10049288061336,
            135.19934282584884,
            -74.00054764512596,
            -184.79956188389923,
            328.,
            -208.30093099671413,
            263.500273822563,
            53.29983570646221,
            4.400328587075575,
            -32.900054764512596,
            -86.69961664841183,
            -124.09994523548741,
            -16.99945235487404,
            3.7995618838992335,
            62.900054764512596,
            60.900054764512596,
            -8.800109529025193,
            -11.,
            83.29983570646222,
            -217.39978094194961,
            2.1998904709748084,
            -57.900054764512596,
            -34.70016429353779,
            58.80010952902519,
            35.600219058050385,
            -90.19989047097481,
            -68.80010952902519,
            70.,
            -54.900054764512596,
            -45.099945235487404,
            0.,
            -13.099945235487404,
            33.900054764512596,
            -10.099945235487404,
            -41.,
            -0.9000547645125958,
            -20.900054764512596,
            28.,
            18.,
            -12.,
            6.,
            -22.,
            11.,
            8.,
            8.,
            -4.,
            -14.099945235487404,
            -9.,
            7.,
            1.,
            -13.,
            2.,
            5.,
            -8.900054764512596,
            16.,
            5.,
            -5.,
            8.,
            -18.,
            8.,
            10.,
            -20.,
            1.,
            14.,
            -11.,
            5.,
            12.,
            -3.,
            1.,
            -2.,
            -2.,
            8.,
            2.,
            10.,
            -0.9000547645125958,
            -2.,
            -1.,
            2.,
            -3.,
            -4.,
            2.,
            2.,
            1.,
            -5.,
            2.,
            -2.,
            6.,
            6.,
            -4.,
            4.,
            0.,
            0.,
            -2.,
            2.,
            4.,
            2.,
            0.,
            0.,
            -6.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
        ];

        let end_coeffs: [f64; 195] = [
            -31519.312979189486,
            -2298.,
            5918.1021358159915,
            -692.2916210295729,
            2911.896221248631,
            -1068.4958926615552,
            959.0807776560789,
            1104.209200438116,
            1026.4975355969332,
            -1476.4958926615552,
            -338.09556407447974,
            1250.9027929901424,
            12.294906900328586,
            590.8896495071194,
            510.1070646221249,
            877.1993428258488,
            632.4975355969332,
            197.3986856516977,
            657.9011500547645,
            -71.3986856516977,
            -366.69687842278205,
            -207.30147864184008,
            137.59802847754656,
            -72.00164293537787,
            -186.3986856516977,
            328.,
            -204.9027929901424,
            262.5008214676889,
            53.89950711938664,
            3.2009857612267254,
            -32.70016429353779,
            -88.09884994523549,
            -124.29983570646222,
            -18.998357064622123,
            5.3986856516977,
            62.70016429353779,
            60.70016429353779,
            -8.400328587075576,
            -11.,
            83.89950711938664,
            -218.19934282584884,
            2.599671412924425,
            -57.70016429353779,
            -34.10049288061336,
            58.40032858707558,
            34.80065717415115,
            -90.59967141292442,
            -68.40032858707558,
            70.,
            -54.70016429353779,
            -45.29983570646221,
            0.,
            -13.299835706462213,
            33.70016429353779,
            -10.299835706462213,
            -41.,
            -0.7001642935377875,
            -20.70016429353779,
            28.,
            18.,
            -12.,
            6.,
            -22.,
            11.,
            8.,
            8.,
            -4.,
            -14.299835706462213,
            -9.,
            7.,
            1.,
            -13.,
            2.,
            5.,
            -8.700164293537787,
            16.,
            5.,
            -5.,
            8.,
            -18.,
            8.,
            10.,
            -20.,
            1.,
            14.,
            -11.,
            5.,
            12.,
            -3.,
            1.,
            -2.,
            -2.,
            8.,
            2.,
            10.,
            -0.7001642935377875,
            -2.,
            -1.,
            2.,
            -3.,
            -4.,
            2.,
            2.,
            1.,
            -5.,
            2.,
            -2.,
            6.,
            6.,
            -4.,
            4.,
            0.,
            0.,
            -2.,
            2.,
            4.,
            2.,
            0.,
            0.,
            -6.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
        ];
        let (a, b) = shval3(59.9, -109.9, 1.1, 10, &start_coeffs, &end_coeffs);
        assert_float_eq!(a.north, 7074.026894642207, rel <= 1e-12);
        assert_float_eq!(a.east, 4596.188334222297, rel <= 1e-12);
        assert_float_eq!(a.down, 62650.49087477549, rel <= 1e-12);
        assert_float_eq!(b.north, 7080.6056093223815, rel <= 1e-12);
        assert_float_eq!(b.east, 4593.347978266618, rel <= 1e-12);
        assert_float_eq!(b.down, 62630.350967492595, rel <= 1e-12);
    }
}

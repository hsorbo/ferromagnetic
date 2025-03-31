use ferromagnetic::{igrf, MagneticComponents, OrthogonalStrength};
use igrf::IGRFresults;
use std::path::Path;

fn lines_to_igrf(chunks: &[&str]) -> IGRFresults {
    IGRFresults {
        result: MagneticComponents {
            declination: chunks[1].parse::<f64>().unwrap(),
            inclination: chunks[3].parse::<f64>().unwrap(),
            horizontal_intensity: chunks[5].parse::<f64>().unwrap(),
            orthogonal_strength: OrthogonalStrength {
                north: chunks[7].parse::<f64>().unwrap(),
                east: chunks[9].parse::<f64>().unwrap(),
                down: chunks[11].parse::<f64>().unwrap(),
            },
            total_intensity: chunks[13].parse::<f64>().unwrap(),
        },
        sv: MagneticComponents {
            declination: chunks[2].parse::<f64>().unwrap(),
            inclination: chunks[4].parse::<f64>().unwrap(),
            horizontal_intensity: chunks[6].parse::<f64>().unwrap(),
            orthogonal_strength: OrthogonalStrength {
                north: chunks[8].parse::<f64>().unwrap(),
                east: chunks[10].parse::<f64>().unwrap(),
                down: chunks[12].parse::<f64>().unwrap(),
            },
            total_intensity: chunks[14].parse::<f64>().unwrap(),
        },
    }
}

const NEAR_POLE_MAX_REL_TOL: f64 = 0.03;
const MAX_REL_TOL: f64 = 0.005;
const NEAR_POLE_TOLERANCE: f64 = 0.001;
const D_I_ABS_TOL: f64 = 0.005;

fn get_max_allowed_relative_error(lat: f64) -> f64 {
    fn is_near_pole(lat: f64) -> bool {
        90.0 - lat.abs() <= NEAR_POLE_TOLERANCE
    }
    if is_near_pole(lat) {
        NEAR_POLE_MAX_REL_TOL
    } else {
        MAX_REL_TOL
    }
}

fn get_max_allowed_absolute_tolerance(value: f64) -> f64 {
    let abs_tol = 0.15;
    // For some values (H, X, Y, Z) that are really small it's hard to calculate the accuracy due to the fact that
    // results from `FORTRAN` are rounded. That's why for some small values (close to the magnetic pole) absolute and relative accuracies
    // are increased.
    if value < 70. {
        abs_tol * 4.
    } else {
        abs_tol
    }
}
fn assert_close(a: f64, b: f64, rel_tol: f64, abs_tol: f64) {
    let abs_diff = (a - b).abs();
    let a_abs = a.abs();
    let b_abs = b.abs();
    let lhs = (rel_tol * a_abs.max(b_abs)).max(abs_tol);
    assert!(
        abs_diff <= lhs,
        "a={} b={} diff={} lhs={}, rel_tol={}, abs_tol={}",
        a,
        b,
        abs_diff,
        lhs,
        rel_tol,
        abs_tol
    );
}

#[test]
fn test_igrf_data() {
    let testdata = Path::new(env!("CARGO_MANIFEST_DIR")).join("testdata/igrf");
    let igrf = igrf::IGRF::default();
    for x in 1..11 {
        let contents = std::fs::read_to_string(testdata.join(format!("set{}", x))).unwrap();
        let raw_hdr = contents
            .lines()
            .next()
            .unwrap()
            .split_whitespace()
            .collect::<Vec<_>>();
        let lat = raw_hdr[1].parse::<f64>().unwrap();
        let lon = raw_hdr[4].parse::<f64>().unwrap();
        let alt = raw_hdr[5].parse::<f64>().unwrap();
        for line in contents.lines().skip(2) {
            let chunks = line.split_whitespace().collect::<Vec<_>>();
            let date = chunks[0].parse::<f64>().unwrap();
            let expected = lines_to_igrf(&chunks).result;
            let actual = igrf.calc(lat, lon, alt, date).result;

            assert_close(
                actual.declination,
                expected.declination,
                get_max_allowed_relative_error(lat),
                D_I_ABS_TOL,
            );

            assert_close(
                actual.inclination,
                expected.inclination,
                get_max_allowed_relative_error(lat),
                D_I_ABS_TOL,
            );

            assert_close(
                actual.horizontal_intensity,
                expected.horizontal_intensity,
                get_max_allowed_relative_error(lat),
                get_max_allowed_absolute_tolerance(actual.horizontal_intensity),
            );

            assert_close(
                actual.orthogonal_strength.north,
                expected.orthogonal_strength.north,
                get_max_allowed_relative_error(lat),
                get_max_allowed_absolute_tolerance(actual.orthogonal_strength.north),
            );

            assert_close(
                actual.orthogonal_strength.east,
                expected.orthogonal_strength.east,
                get_max_allowed_relative_error(lat),
                get_max_allowed_absolute_tolerance(actual.orthogonal_strength.east),
            );

            assert_close(
                actual.orthogonal_strength.down,
                expected.orthogonal_strength.down,
                get_max_allowed_relative_error(lat),
                get_max_allowed_absolute_tolerance(actual.orthogonal_strength.down),
            );

            assert_close(
                actual.total_intensity,
                expected.total_intensity,
                get_max_allowed_relative_error(lat),
                get_max_allowed_absolute_tolerance(actual.total_intensity),
            );
        }
    }
}

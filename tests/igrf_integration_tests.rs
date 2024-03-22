use ferromagnetic::igrf;
use igrf::IGRFresults;
use std::path::Path;

fn lines_to_igrf(chunks: &[&str]) -> IGRFresults {
    IGRFresults {
        declination: chunks[1].parse::<f64>().unwrap(),
        declination_sv: chunks[2].parse::<f64>().unwrap(),
        inclination: chunks[3].parse::<f64>().unwrap(),
        inclination_sv: chunks[4].parse::<f64>().unwrap(),
        horizontal_intensity: chunks[5].parse::<f64>().unwrap(),
        horizontal_sv: chunks[6].parse::<f64>().unwrap(),
        north_component: chunks[7].parse::<f64>().unwrap(),
        north_sv: chunks[8].parse::<f64>().unwrap(),
        east_component: chunks[9].parse::<f64>().unwrap(),
        east_sv: chunks[10].parse::<f64>().unwrap(),
        vertical_component: chunks[11].parse::<f64>().unwrap(),
        vertical_sv: chunks[12].parse::<f64>().unwrap(),
        total_intensity: chunks[13].parse::<f64>().unwrap(),
        total_sv: chunks[14].parse::<f64>().unwrap(),
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
            let expected = lines_to_igrf(&chunks);
            let actual = igrf.calc(lat, lon, alt, date);

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
                actual.north_component,
                expected.north_component,
                get_max_allowed_relative_error(lat),
                get_max_allowed_absolute_tolerance(actual.north_component),
            );

            assert_close(
                actual.east_component,
                expected.east_component,
                get_max_allowed_relative_error(lat),
                get_max_allowed_absolute_tolerance(actual.east_component),
            );

            assert_close(
                actual.vertical_component,
                expected.vertical_component,
                get_max_allowed_relative_error(lat),
                get_max_allowed_absolute_tolerance(actual.vertical_component),
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

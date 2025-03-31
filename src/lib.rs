pub mod igrf;

pub struct OrthogonalStrength {
    /// North component (X) (nT)
    pub north: f64,
     /// East component (Y) (nT)
    pub east: f64,
    /// Down / Vertical component (Z) (nT)
    pub down: f64,
}

impl Default for OrthogonalStrength {
    fn default() -> Self {
        OrthogonalStrength {
            north: 0.0,
            east: 0.0,
            down: 0.0,
        }
    }
}

/// Components needed to determine Earth's magnetic field at a given location
/// Detailed info: https://www.geomag.nrcan.gc.ca/mag_fld/comp-en.php
pub struct MagneticComponents {
    /// Declination (D) (degrees)
    pub declination: f64,
    /// Inclination (D) (degrees)
    pub inclination: f64,
    /// Horizontal intensity (H) (nT)
    pub horizontal_intensity: f64,
    /// Orthogonal strength components
    pub orthogonal_strength: OrthogonalStrength,
    /// Total intensity (F) (nT)
    pub total_intensity: f64,
}
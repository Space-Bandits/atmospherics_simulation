use crate::mixture::Mixture;

pub struct FluidVolume {
    mixture: Mixture,
    /// Litres
    volume: f32,
}

pub struct VolumeProperties {
    /// Litres
    liquid_volume: f32,
    /// Kelvin
    temperature: f32,
    /// kPa
    pressure: f32,
}

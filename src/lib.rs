pub mod flow;
pub mod fluid_properties;
pub mod fluid_volume;
pub mod mixture;

pub const IDEAL_GAS_CONSTANT: f32 = 8.3144;

pub trait ToKelvin {
    fn to_kelvin(self) -> Self;
}

pub trait ToCelcius {
    fn to_celcius(self) -> Self;
}

impl ToKelvin for f32 {
    fn to_kelvin(self) -> Self {
        self + 273.15
    }
}

impl ToCelcius for f32 {
    fn to_celcius(self) -> Self {
        self - 273.15
    }
}

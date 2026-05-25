use thiserror::Error;

/// Properties for a type of fluid.
pub struct FluidProperties {
    pub fluid_type: FluidPropertiesType,
    /// Joules per mol per kelvin
    pub heat_capactity: f32,
}

pub enum FluidPropertiesType {
    Liquid(LiquidFluidProperties),
    Gas(GasFluidProperties),
}

/// Properties that a liquid has.
///
/// belongs to a [FluidProperties].
pub struct LiquidFluidProperties {
    /// Litres per mol
    pub density: f32,
}

/// Properties that a gas has.
///
/// belongs to a [FluidProperties].
pub struct GasFluidProperties {}

/// A collection of [FluidProperties] and their unique ids.
#[derive(Default)]
pub struct FluidCollection {
    fluids: Vec<FluidProperties>,
}

impl FluidCollection {
    /// Add a new fluid and assign it an id.
    pub fn add_fluid(&mut self, properties: FluidProperties) -> usize {
        let fluid_id = self.fluids.len();

        self.fluids.push(properties);

        fluid_id
    }

    /// Get the properties of a particular fluid
    pub fn get_properties(&self, fluid_id: usize) -> Result<&FluidProperties, InvalidFluidId> {
        self.fluids
            .get(fluid_id)
            .ok_or(InvalidFluidId(fluid_id, self.fluids.len()))
    }

    /// Returns an iterator over fluid properties and their ids.
    pub fn iter(&self) -> impl Iterator<Item = (usize, &FluidProperties)> {
        self.fluids.iter().enumerate()
    }
}

#[derive(Error, Debug)]
#[error(
    "This fluid id index is out of range, asked for fluid index {0} but there are only {1} fluids"
)]
pub struct InvalidFluidId(usize, usize);

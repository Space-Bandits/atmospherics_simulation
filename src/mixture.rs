use crate::fluid_properties::{FluidCollection, FluidTypeProperties, InvalidFluidId};

/// An amount of a specific fluid.
#[derive(Clone, Copy, Debug)]
pub struct Fluid {
    pub fluid_id: usize,
    /// Moles
    pub moles: f32,
}

/// A mixture of fluids with a temperature.
#[derive(Default)]
pub struct Mixture {
    fluids: Vec<Fluid>,
    /// Joules
    energy: f32,
}

/// The properties of a [Mixture] that are derived from its makeup.
#[derive(Clone, Copy, Debug)]
pub struct ComputedMixtureProperties {
    /// Joules per Kelvin
    pub heat_capacity: f32,
    /// Litres
    pub liquid_volume: f32,
    /// Moles of ideal gas
    pub gas_moles: f32,
    /// Kelvin
    pub temperature: f32,
}

impl Mixture {
    pub fn add_fluid(&mut self, fluid: Fluid) {
        match self
            .fluids
            .binary_search_by(|probe| probe.fluid_id.cmp(&fluid.fluid_id))
        {
            Ok(index) => self.fluids[index].moles += fluid.moles,
            Err(index) => self.fluids.insert(index, fluid),
        }
    }

    pub fn add_fluid_at_temperature(
        &mut self,
        collection: &FluidCollection,
        fluid: Fluid,
        temperature: f32,
    ) -> Result<(), InvalidFluidId> {
        let properties = collection.get_properties(fluid.fluid_id)?;

        self.energy += fluid.moles * properties.heat_capactity * temperature;

        self.add_fluid(fluid);

        Ok(())
    }

    pub fn from_fluid_at_temperature(
        collection: &FluidCollection,
        fluid: Fluid,
        temperature: f32,
    ) -> Result<Self, InvalidFluidId> {
        let mut mixture = Mixture {
            fluids: Vec::with_capacity(1),
            energy: 0.,
        };

        mixture.add_fluid_at_temperature(collection, fluid, temperature)?;

        Ok(mixture)
    }

    pub fn add_mixture(&mut self, mixture: Self) {
        for fluid in mixture.fluids {
            self.add_fluid(fluid);
        }
    }

    pub fn add_energy(&mut self, energy: f32) {
        self.energy += energy;
    }

    pub fn energy(&self) -> f32 {
        self.energy
    }

    pub fn total_moles(&self) -> f32 {
        self.fluids.iter().map(|fluid| fluid.moles).sum()
    }

    pub fn compute_mixture_properties(
        &self,
        collection: &FluidCollection,
    ) -> Result<ComputedMixtureProperties, InvalidFluidId> {
        let mut heat_capacity = 0.;
        let mut liquid_volume = 0.;
        let mut gas_moles = 0.;

        for fluid in &self.fluids {
            let properties = collection.get_properties(fluid.fluid_id)?;

            heat_capacity += properties.heat_capactity * fluid.moles;

            match properties.fluid_type {
                FluidTypeProperties::Liquid(ref liquid_properties) => {
                    liquid_volume += liquid_properties.density * fluid.moles;
                }
                FluidTypeProperties::Gas(_) => {
                    gas_moles += fluid.moles;
                }
            }
        }

        let temperature = heat_capacity * self.energy;

        Ok(ComputedMixtureProperties {
            heat_capacity,
            liquid_volume,
            gas_moles,
            temperature,
        })
    }
}

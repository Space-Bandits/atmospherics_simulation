use crate::fluid_properties::{
    FluidCollection, FluidProperties, FluidTypeProperties, InvalidFluidId,
};

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
#[derive(Clone, Copy, Debug, Default)]
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
    pub fn new(fluids: impl IntoIterator<Item = Fluid>, energy: f32) -> Self {
        let mut mixture = Mixture {
            fluids: Vec::new(),
            energy,
        };

        for fluid in fluids {
            mixture.add_fluid(fluid);
        }

        mixture
    }

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

        self.energy += fluid.moles * properties.molar_heat_capactity * temperature;

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

    pub fn from_fluids_at_temperature(
        collection: &FluidCollection,
        fluids: impl IntoIterator<Item = Fluid>,
        temperature: f32,
    ) -> Result<Self, InvalidFluidId> {
        let mut mixture = Mixture {
            fluids: Vec::new(),
            energy: 0.,
        };

        for fluid in fluids {
            mixture.add_fluid_at_temperature(collection, fluid, temperature)?;
        }

        Ok(mixture)
    }

    pub fn add_mixture(&mut self, mixture: Self) {
        // Merge the two sorted lists of fluids together, combine duplicate entries.

        let mut probe_index = 0;

        for fluid in mixture.fluids {
            loop {
                if let Some(probed_fluid) = self.fluids.get_mut(probe_index) {
                    match fluid.fluid_id.cmp(&probed_fluid.fluid_id) {
                        std::cmp::Ordering::Less => {
                            // If the id of the added fluid is before the probed fluid insert and shift up at the probe.
                            self.fluids.insert(probe_index, fluid);
                            // Shift probe up after insert
                            probe_index += 1;
                            break;
                        }
                        std::cmp::Ordering::Equal => {
                            probed_fluid.moles += fluid.moles;
                            // There won't be two of the same fluid so shift probe up.
                            probe_index += 1;
                            break;
                        }
                        std::cmp::Ordering::Greater => {
                            // This fluid should be inserted after the probed fluid, shift the probe and check next.
                            probe_index += 1;
                            continue;
                        }
                    }
                }

                self.fluids.push(fluid);
                // Shift probe up after insert
                probe_index += 1;
                break;
            }
        }

        self.energy += mixture.energy;
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

            heat_capacity += properties.molar_heat_capactity * fluid.moles;

            match properties.fluid_type {
                FluidTypeProperties::Liquid(ref liquid_properties) => {
                    liquid_volume += liquid_properties.density * fluid.moles;
                }
                FluidTypeProperties::Gas(_) => {
                    gas_moles += fluid.moles;
                }
            }
        }

        let temperature = self.energy / heat_capacity;

        Ok(ComputedMixtureProperties {
            heat_capacity,
            liquid_volume,
            gas_moles,
            temperature,
        })
    }

    /// Selectively extracts fluids and energy from the mixture.
    ///
    /// If an error is returned the mixture may still have had some fluids removed.
    pub fn extract_fluids(
        &mut self,
        collection: &FluidCollection,
        mut predicate: impl FnMut(&Fluid, &FluidProperties) -> f32,
    ) -> Result<Self, InvalidFluidId> {
        let mut mixture = Mixture::default();

        // Joules per Kelvin
        let mut heat_capacity_kept = 0.;
        let mut heat_capacity_extracted = 0.;

        for fluid in &mut self.fluids {
            let fluid_properties = collection.get_properties(fluid.fluid_id)?;

            let extract_moles = predicate(fluid, fluid_properties).min(fluid.moles);

            if extract_moles <= 0. {
                continue;
            }

            fluid.moles -= extract_moles;

            // Ordering is retained
            mixture.fluids.push(Fluid {
                fluid_id: fluid.fluid_id,
                moles: extract_moles,
            });

            heat_capacity_kept += fluid_properties.molar_heat_capactity * fluid.moles;
            heat_capacity_extracted += fluid_properties.molar_heat_capactity * extract_moles;
        }

        let heat_capacity = heat_capacity_kept + heat_capacity_extracted;

        mixture.energy = self.energy * (heat_capacity_extracted / heat_capacity);
        self.energy = self.energy * (heat_capacity_kept / heat_capacity);

        Ok(mixture)
    }
}

#[cfg(test)]
mod tests {
    use crate::mixture::{Fluid, Mixture};

    #[test]
    fn add_mixture_identical_fluid() {
        let mut mixture = Mixture {
            fluids: vec![Fluid {
                fluid_id: 5,
                moles: 1.,
            }],
            energy: 10.,
        };

        mixture.add_mixture(Mixture {
            fluids: vec![Fluid {
                fluid_id: 5,
                moles: 1.,
            }],
            energy: 10.,
        });

        assert_eq!(mixture.fluids[0].moles, 1. + 1.);
    }

    #[test]
    fn add_mixture_lower_fluid() {
        let mut mixture = Mixture {
            fluids: vec![Fluid {
                fluid_id: 5,
                moles: 1.,
            }],
            energy: 10.,
        };

        mixture.add_mixture(Mixture {
            fluids: vec![Fluid {
                fluid_id: 4,
                moles: 1.,
            }],
            energy: 10.,
        });

        assert_eq!(mixture.fluids[0].fluid_id, 4);
        assert_eq!(mixture.fluids[1].fluid_id, 5);
    }

    #[test]
    fn add_mixture_higher_fluid() {
        let mut mixture = Mixture {
            fluids: vec![Fluid {
                fluid_id: 5,
                moles: 1.,
            }],
            energy: 10.,
        };

        mixture.add_mixture(Mixture {
            fluids: vec![Fluid {
                fluid_id: 6,
                moles: 1.,
            }],
            energy: 10.,
        });

        assert_eq!(mixture.fluids[0].fluid_id, 5);
        assert_eq!(mixture.fluids[1].fluid_id, 6);
    }
}

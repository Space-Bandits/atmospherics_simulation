use crate::{
    IDEAL_GAS_CONSTANT,
    fluid_properties::{FluidCollection, InvalidFluidId},
    mixture::{ComputedMixtureProperties, Mixture},
};

pub struct FluidVolume {
    mixture: Mixture,
    /// Litres
    volume: f32,
}

/// The properties of a [FluidVolume] that are derived from its makeup.
#[derive(Clone, Copy, Debug)]
pub struct ComputedVolumeProperties {
    pub mixture_properties: ComputedMixtureProperties,
    /// kPa
    pub pressure: f32,
}

impl FluidVolume {
    pub fn new(volume: f32, mixture: Mixture) -> Self {
        FluidVolume { mixture, volume }
    }

    pub fn new_empty(volume: f32) -> Self {
        FluidVolume {
            mixture: Mixture::default(),
            volume,
        }
    }

    pub fn calculate_properties(
        &self,
        collection: &FluidCollection,
    ) -> Result<ComputedVolumeProperties, InvalidFluidId> {
        let mixture_properties = self.mixture.compute_mixture_properties(collection)?;

        let gas_volume = self.volume - mixture_properties.liquid_volume;

        let pressure =
            (mixture_properties.gas_moles * IDEAL_GAS_CONSTANT * mixture_properties.temperature)
                / gas_volume;

        Ok(ComputedVolumeProperties {
            mixture_properties,
            pressure,
        })
    }
}

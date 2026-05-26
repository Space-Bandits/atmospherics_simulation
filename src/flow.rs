use std::collections::HashMap;

use thiserror::Error;

use crate::{
    flow::FlowSimulationError::InvalidVolumeId,
    fluid_properties::{FluidCollection, InvalidFluidId},
    fluid_volume::{ComputedVolumeProperties, FluidVolume},
    mixture::ComputedMixtureProperties,
};

pub struct FlowSimulationState<K> {
    volumes: Vec<K>,
    volume_states: HashMap<K, FlowVolumeNodeSimulationState>,
}

struct FlowVolumeNodeSimulationState {
    computed_properties: ComputedVolumeProperties,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum EdgeFluidPriority {
    Equal,
    FavorLiquid,
    FavorGas,
}

/// Trait methods that allow a [FlowSimulationState] to query and modify [FluidVolume] states.
///
/// - `K` is the key that indentifies a volume.
/// - `E` is an error type that gets bubbled up through [FlowSimulationState::step] if any of the methods fail.
pub trait SimulationQueries<K, E> {
    /// Should return an iterator of all other edges that fluid can flow into from a volume.
    ///
    /// - The `K` is the id of the connected volume.
    /// - The `f32` is the cross sectional area of the connection. A larger value means more flow.
    /// - The `EdgeFluidPriority` determines whether fluid, gas or both will be pushed through the connection first.
    ///
    /// Should always return the exact same edges in the same order for any volume within one simulation step.
    fn get_edges(
        &self,
        volume_id: &K,
    ) -> Result<impl IntoIterator<Item = (K, f32, EdgeFluidPriority)>, E>;

    fn get_volume(&self, volume_id: &K) -> Result<&FluidVolume, E>;

    fn get_volume_mut(&mut self, volume_id: &K) -> Result<&mut FluidVolume, E>;
}

impl<K> Default for FlowSimulationState<K> {
    fn default() -> Self {
        FlowSimulationState {
            volumes: Vec::new(),
            volume_states: HashMap::new(),
        }
    }
}

impl<K> FlowSimulationState<K>
where
    K: std::hash::Hash + Eq + Copy,
{
    pub fn add_volume(&mut self, volume_id: K) {
        self.volumes.push(volume_id);

        self.volume_states.insert(
            volume_id,
            FlowVolumeNodeSimulationState {
                computed_properties: ComputedVolumeProperties {
                    mixture_properties: ComputedMixtureProperties {
                        heat_capacity: 0.,
                        liquid_volume: 0.,
                        gas_moles: 0.,
                        temperature: 0.,
                    },
                    pressure: 0.,
                },
            },
        );
    }

    pub fn step<E>(
        &mut self,
        collection: &FluidCollection,
        mut queries: impl SimulationQueries<K, E>,
    ) -> Result<(), FlowSimulationError<K, E>> {
        for volume_id in &self.volumes {
            self.volume_states
                .get_mut(volume_id)
                .expect("This volume id should exist")
                .computed_properties = queries
                .get_volume(volume_id)
                .map_err(|err| FlowQueryError(err))?
                .calculate_properties(collection)?
        }

        for volume_id in &self.volumes {
            let volume_properties = self
                .volume_states
                .get(volume_id)
                .expect("This volume id should exist")
                .computed_properties;

            for (other_volume_id, cross_section, priority) in queries
                .get_edges(volume_id)
                .map_err(|err| FlowQueryError(err))?
            {
                let other_volume_properties = self
                    .volume_states
                    .get(&other_volume_id)
                    .ok_or(FlowSimulationError::InvalidVolumeId(other_volume_id))?
                    .computed_properties;
            }
        }

        Ok(())
    }
}

#[derive(Error, Debug)]
#[error("Encountered data query error whilst stepping fluid simulation: {0}")]
pub struct FlowQueryError<E>(pub E);

#[derive(Error, Debug)]
pub enum FlowSimulationError<K, E> {
    #[error("{0}")]
    QueryError(#[from] FlowQueryError<E>),
    #[error("Encountered invalid fluid id whilst stepping fluid simulation: {0}")]
    InvalidFluidId(#[from] InvalidFluidId),
    #[error("Volume id did not exist in the simulation: {0}")]
    InvalidVolumeId(K),
}

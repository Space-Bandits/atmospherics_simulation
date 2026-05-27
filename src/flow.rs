use std::collections::HashMap;

use thiserror::Error;

use crate::{
    fluid_properties::{FluidCollection, FluidTypeProperties, InvalidFluidId},
    fluid_volume::{ComputedVolumeProperties, FluidVolume},
    mixture::Mixture,
};

pub struct FlowSimulationState<K> {
    volumes: Vec<K>,
    volume_states: HashMap<K, FlowVolumeNodeSimulationState>,
    edges: Vec<FlowEdgeState>,
}

#[derive(Default)]
struct FlowVolumeNodeSimulationState {
    volume: f32,
    computed_properties: ComputedVolumeProperties,
    /// Indices of edges sending liquid from this volume.
    send_edges: Vec<usize>,
    /// The maximum volume of liquid that an edge would remove from this volume independently.
    send_max_volume: f32,
    /// Indices of edges sending liquid from to volume.
    receive_edges: Vec<usize>,
    /// The maximum volume of liquid that an edge would add to this volume independently.
    receive_max_volume: f32,
}

struct FlowEdgeState {
    move_liquid: f32,
    limited_move_liquid: f32,
    extracted_mixture: Mixture,
}

/// Trait methods that allow a [FlowSimulationState] to query and modify [FluidVolume] states.
///
/// - `K` is the key that indentifies a volume.
/// - `E` is an error type that gets bubbled up through [FlowSimulationState::step] if any of the methods fail.
pub trait SimulationQueries<K, E> {
    /// Should return an iterator of all other edges that liquid can flow into from a volume.
    ///
    /// Should always return the exact same edges in the same order for any volume within one simulation step.
    fn get_liquid_edges(&self, volume_id: &K) -> Result<impl IntoIterator<Item = K>, E>;

    /// Should return an iterator of all other edges that gas can flow into from a volume.
    ///
    /// Should always return the exact same edges in the same order for any volume within one simulation step.
    fn get_gas_edges(&self, volume_id: &K) -> Result<impl IntoIterator<Item = K>, E>;

    fn get_volume(&self, volume_id: &K) -> Result<&FluidVolume, E>;

    fn get_volume_mut(&mut self, volume_id: &K) -> Result<&mut FluidVolume, E>;
}

impl<K> Default for FlowSimulationState<K> {
    fn default() -> Self {
        FlowSimulationState {
            volumes: Vec::new(),
            volume_states: HashMap::new(),
            edges: Vec::new(),
        }
    }
}

impl<K> FlowSimulationState<K>
where
    K: std::hash::Hash + Eq + Copy,
{
    pub fn add_volume(&mut self, volume_id: K) {
        self.volumes.push(volume_id);

        self.volume_states
            .insert(volume_id, FlowVolumeNodeSimulationState::default());
    }

    pub fn step<E>(
        &mut self,
        collection: &FluidCollection,
        queries: &mut impl SimulationQueries<K, E>,
    ) -> Result<(), FlowSimulationError<K, E>> {
        self.step_liquid(collection, queries)?;
        self.step_gas(collection, queries)?;

        Ok(())
    }

    fn step_liquid<E>(
        &mut self,
        collection: &FluidCollection,
        queries: &mut impl SimulationQueries<K, E>,
    ) -> Result<(), FlowSimulationError<K, E>> {
        self.edges.clear();

        for volume_id in &self.volumes {
            let volume = queries
                .get_volume(volume_id)
                .map_err(|err| FlowQueryError(err))?;

            let state = self.volume_states.get_mut(volume_id).unwrap();

            *state = FlowVolumeNodeSimulationState {
                volume: volume.volume(),
                computed_properties: volume.calculate_properties(collection)?,
                send_edges: Vec::new(),
                send_max_volume: 0.,
                receive_edges: Vec::new(),
                receive_max_volume: 0.,
            };
        }

        // Pass over all edges and determine what amount of liquid each of them would move in isolation
        for volume_id in &self.volumes {
            let &FlowVolumeNodeSimulationState {
                volume,
                computed_properties,
                ..
            } = self.volume_states.get(volume_id).unwrap();

            for other_volume_id in queries
                .get_liquid_edges(volume_id)
                .map_err(|err| FlowQueryError(err))?
            {
                let &FlowVolumeNodeSimulationState {
                    volume: other_volume,
                    computed_properties: other_computed_properties,
                    ..
                } = self
                    .volume_states
                    .get(&other_volume_id)
                    .ok_or(FlowSimulationError::InvalidVolumeId(other_volume_id))?;

                // What volume of liquid needs to be moved into the other tank in order to equalize the fill levels
                let move_liquid = {
                    let v1 = volume;
                    let v2 = other_volume;
                    let l1 = computed_properties.mixture_properties.liquid_volume;
                    let l2 = other_computed_properties.mixture_properties.liquid_volume;

                    (l1 * v2 - l2 * v1) / (v1 + v2)
                };

                if move_liquid < 0. {
                    continue;
                }

                let edge_index = self.edges.len();

                self.edges.push(FlowEdgeState {
                    move_liquid: move_liquid,
                    limited_move_liquid: move_liquid,
                    extracted_mixture: Mixture::default(),
                });

                let state = self.volume_states.get_mut(volume_id).unwrap();
                state.send_edges.push(edge_index);
                state.send_max_volume = state.send_max_volume.max(move_liquid);

                let state = self.volume_states.get_mut(&other_volume_id).unwrap();
                state.receive_edges.push(edge_index);
                state.receive_max_volume = state.receive_max_volume.max(move_liquid);
            }
        }

        // Flow is limited using this process for both outgoing and incomming liquid:
        //
        // Find the sum of all liquid flow into/out of a volume, as well as the maximum inward and outward flow.
        // Scale all edge flows uniformly so that the flow into/out of the volume is the maximum.
        //
        // This stops situations where multiple edges compounding on each other try to overfill or overdrain a volume because they don't consider the whole graph.

        for volume_id in &self.volumes {
            let state = self.volume_states.get(volume_id).unwrap();

            let total_send_volume: f32 = state
                .send_edges
                .iter()
                .map(|&edge_index| self.edges.get(edge_index).unwrap().move_liquid)
                .sum();

            let scale = state.send_max_volume / total_send_volume;

            if scale.is_finite() {
                for &edge_index in &state.send_edges {
                    let edge = self.edges.get_mut(edge_index).unwrap();

                    edge.limited_move_liquid = edge.move_liquid * scale;
                }
            }

            let total_receive_volume: f32 = state
                .receive_edges
                .iter()
                .map(|&edge_index| self.edges.get(edge_index).unwrap().move_liquid)
                .sum();

            let scale = state.receive_max_volume / total_receive_volume;

            if scale.is_finite() {
                for &edge_index in &state.receive_edges {
                    let edge = self.edges.get_mut(edge_index).unwrap();

                    edge.limited_move_liquid =
                        edge.limited_move_liquid.min(edge.move_liquid * scale);
                }
            }
        }

        // Extract fluids from mixtures

        for volume_id in &self.volumes {
            let state = self.volume_states.get(volume_id).unwrap();

            let volume = queries
                .get_volume_mut(volume_id)
                .map_err(|err| FlowQueryError(err))?;

            for &edge_index in &state.send_edges {
                let edge = self.edges.get_mut(edge_index).unwrap();

                edge.extracted_mixture =
                    volume.mixture.extract_fluids(collection, |_, properties| {
                        if let FluidTypeProperties::Liquid(liquid_properties) =
                            &properties.fluid_type
                        {
                            edge.limited_move_liquid / liquid_properties.density
                        } else {
                            0.
                        }
                    })?;
            }
        }

        // Insert mixtures

        for volume_id in &self.volumes {
            let state = self.volume_states.get(volume_id).unwrap();

            let volume = queries
                .get_volume_mut(volume_id)
                .map_err(|err| FlowQueryError(err))?;

            for &edge_index in &state.receive_edges {
                let edge = self.edges.get_mut(edge_index).unwrap();

                volume
                    .mixture
                    .add_mixture(std::mem::take(&mut edge.extracted_mixture));
            }
        }

        Ok(())
    }

    fn step_gas<E>(
        &mut self,
        collection: &FluidCollection,
        queries: &mut impl SimulationQueries<K, E>,
    ) -> Result<(), FlowSimulationError<K, E>> {
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

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
    liquid_edges: Vec<LiquidEdgeState>,
    gas_edges: Vec<GasEdgeState>,
}

#[derive(Default)]
struct FlowVolumeNodeSimulationState {
    volume: f32,
    computed_properties: ComputedVolumeProperties,

    /// Indices of edges sending liquid from this volume.
    send_liquid_edges: Vec<usize>,
    /// Indices of edges sending liquid from to volume.
    receive_liquid_edges: Vec<usize>,
    /// The maximum volume of liquid that an edge would remove from this volume independently.
    send_max_volume: f32,
    /// The maximum volume of liquid that an edge would add to this volume independently.
    receive_max_volume: f32,

    /// Indices of edges sending gas from this volume.
    send_gas_edges: Vec<usize>,
    /// Indices of edges sending gas from to volume.
    receive_gas_edges: Vec<usize>,
    /// The maximum pressure volumes of gas that an edge would remove from this volume independently.
    send_max_pv: f32,
    /// The maximum pressure volumes of gas that an edge would add to this volume independently.
    receive_max_pv: f32,
}

struct LiquidEdgeState {
    move_volume: f32,
    limited_move_volume: f32,
    extracted_mixture: Mixture,
}

struct GasEdgeState {
    move_pv: f32,
    limited_move_pv: f32,
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
            liquid_edges: Vec::new(),
            gas_edges: Vec::new(),
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
        // Flow is limited using this process for both outgoing and incomming liquid and gas:
        //
        // Find the sum of all fluid flow into/out of a volume, as well as the maximum inward and outward flow.
        // Scale all edge flows uniformly so that the flow into/out of the volume is the maximum.
        //
        // This stops situations where multiple edges compounding on each other try to overfill or overdrain a volume because they don't consider the whole graph.
        //
        // Liquid flow is calculated first and displaces any gas in its way, which is calculated at a second priority after.
        // This means that gas cannot be pushed around by liquids

        self.liquid_edges.clear();
        self.gas_edges.clear();

        // Get volumes, calculate derived properties and clear working values.
        for volume_id in &self.volumes {
            let volume = queries
                .get_volume(volume_id)
                .map_err(|err| FlowQueryError(err))?;

            let state = self.volume_states.get_mut(volume_id).unwrap();

            state.volume = volume.volume();
            state.computed_properties = volume.calculate_properties(collection)?;

            state.send_liquid_edges.clear();
            state.receive_liquid_edges.clear();
            state.send_max_volume = 0.;
            state.receive_max_volume = 0.;

            state.send_gas_edges.clear();
            state.receive_gas_edges.clear();
            state.send_max_pv = 0.;
            state.receive_max_pv = 0.;
        }

        self.step_liquid(collection, queries)?;

        // Recalculate derived values after moving liquids
        for volume_id in &self.volumes {
            let volume = queries
                .get_volume(volume_id)
                .map_err(|err| FlowQueryError(err))?;

            let state = self.volume_states.get_mut(volume_id).unwrap();

            state.computed_properties = volume.calculate_properties(collection)?;
        }

        self.step_gas(collection, queries)?;

        Ok(())
    }

    fn step_liquid<E>(
        &mut self,
        collection: &FluidCollection,
        queries: &mut impl SimulationQueries<K, E>,
    ) -> Result<(), FlowSimulationError<K, E>> {
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
                let move_volume = {
                    let v1 = volume;
                    let v2 = other_volume;
                    let l1 = computed_properties.mixture_properties.liquid_volume;
                    let l2 = other_computed_properties.mixture_properties.liquid_volume;

                    (l1 * v2 - l2 * v1) / (v1 + v2)
                };

                if move_volume < 0. {
                    continue;
                }

                let edge_index = self.liquid_edges.len();

                self.liquid_edges.push(LiquidEdgeState {
                    move_volume,
                    limited_move_volume: move_volume,
                    extracted_mixture: Mixture::default(),
                });

                let state = self.volume_states.get_mut(volume_id).unwrap();
                state.send_liquid_edges.push(edge_index);
                state.send_max_volume = state.send_max_volume.max(move_volume);

                let state = self.volume_states.get_mut(&other_volume_id).unwrap();
                state.receive_liquid_edges.push(edge_index);
                state.receive_max_volume = state.receive_max_volume.max(move_volume);
            }
        }

        for volume_id in &self.volumes {
            let state = self.volume_states.get(volume_id).unwrap();

            let total_send_volume: f32 = state
                .send_liquid_edges
                .iter()
                .map(|&edge_index| self.liquid_edges.get(edge_index).unwrap().move_volume)
                .sum();

            let scale = state.send_max_volume / total_send_volume;

            if scale.is_finite() {
                for &edge_index in &state.send_liquid_edges {
                    let edge = self.liquid_edges.get_mut(edge_index).unwrap();

                    edge.limited_move_volume = edge.move_volume * scale;
                }
            }

            let total_receive_volume: f32 = state
                .receive_liquid_edges
                .iter()
                .map(|&edge_index| self.liquid_edges.get(edge_index).unwrap().move_volume)
                .sum();

            let scale = state.receive_max_volume / total_receive_volume;

            if scale.is_finite() {
                for &edge_index in &state.receive_liquid_edges {
                    let edge = self.liquid_edges.get_mut(edge_index).unwrap();

                    edge.limited_move_volume =
                        edge.limited_move_volume.min(edge.move_volume * scale);
                }
            }
        }

        // Extract fluids from mixtures

        for volume_id in &self.volumes {
            let state = self.volume_states.get(volume_id).unwrap();

            let volume = queries
                .get_volume_mut(volume_id)
                .map_err(|err| FlowQueryError(err))?;

            for &edge_index in &state.send_liquid_edges {
                let edge = self.liquid_edges.get_mut(edge_index).unwrap();

                edge.extracted_mixture =
                    volume.mixture.extract_fluids(collection, |_, properties| {
                        if let FluidTypeProperties::Liquid(liquid_properties) =
                            &properties.fluid_type
                        {
                            edge.limited_move_volume / liquid_properties.density
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

            for &edge_index in &state.receive_liquid_edges {
                let edge = self.liquid_edges.get_mut(edge_index).unwrap();

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
        // Pass over all edges and determine what amount of gas each of them would move in isolation
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
                let move_pv = {
                    let v1 = volume - computed_properties.mixture_properties.liquid_volume;
                    let v2 =
                        other_volume - other_computed_properties.mixture_properties.liquid_volume;
                    let p1 = computed_properties.pressure;
                    let p2 = other_computed_properties.pressure;

                    (p1 * v2 - p2 * v1) / (v1 + v2)
                };

                if move_pv < 0. {
                    continue;
                }

                let edge_index = self.gas_edges.len();

                self.gas_edges.push(GasEdgeState {
                    move_pv,
                    limited_move_pv: move_pv,
                    extracted_mixture: Mixture::default(),
                });

                let state = self.volume_states.get_mut(volume_id).unwrap();
                state.send_gas_edges.push(edge_index);
                state.send_max_volume = state.send_max_volume.max(move_pv);

                let state = self.volume_states.get_mut(&other_volume_id).unwrap();
                state.receive_gas_edges.push(edge_index);
                state.receive_max_volume = state.receive_max_volume.max(move_pv);
            }
        }

        // for volume_id in &self.volumes {
        //     let state = self.volume_states.get(volume_id).unwrap();

        //     let total_send_volume: f32 = state
        //         .send_liquid_edges
        //         .iter()
        //         .map(|&edge_index| self.liquid_edges.get(edge_index).unwrap().move_volume)
        //         .sum();

        //     let scale = state.send_max_volume / total_send_volume;

        //     if scale.is_finite() {
        //         for &edge_index in &state.send_liquid_edges {
        //             let edge = self.liquid_edges.get_mut(edge_index).unwrap();

        //             edge.limited_move_volume = edge.move_volume * scale;
        //         }
        //     }

        //     let total_receive_volume: f32 = state
        //         .receive_liquid_edges
        //         .iter()
        //         .map(|&edge_index| self.liquid_edges.get(edge_index).unwrap().move_volume)
        //         .sum();

        //     let scale = state.receive_max_volume / total_receive_volume;

        //     if scale.is_finite() {
        //         for &edge_index in &state.receive_liquid_edges {
        //             let edge = self.liquid_edges.get_mut(edge_index).unwrap();

        //             edge.limited_move_volume =
        //                 edge.limited_move_volume.min(edge.move_volume * scale);
        //         }
        //     }
        // }

        // // Extract fluids from mixtures

        // for volume_id in &self.volumes {
        //     let state = self.volume_states.get(volume_id).unwrap();

        //     let volume = queries
        //         .get_volume_mut(volume_id)
        //         .map_err(|err| FlowQueryError(err))?;

        //     for &edge_index in &state.send_liquid_edges {
        //         let edge = self.liquid_edges.get_mut(edge_index).unwrap();

        //         edge.extracted_mixture =
        //             volume.mixture.extract_fluids(collection, |_, properties| {
        //                 if let FluidTypeProperties::Liquid(liquid_properties) =
        //                     &properties.fluid_type
        //                 {
        //                     edge.limited_move_volume / liquid_properties.density
        //                 } else {
        //                     0.
        //                 }
        //             })?;
        //     }
        // }

        // // Insert mixtures

        // for volume_id in &self.volumes {
        //     let state = self.volume_states.get(volume_id).unwrap();

        //     let volume = queries
        //         .get_volume_mut(volume_id)
        //         .map_err(|err| FlowQueryError(err))?;

        //     for &edge_index in &state.receive_liquid_edges {
        //         let edge = self.liquid_edges.get_mut(edge_index).unwrap();

        //         volume
        //             .mixture
        //             .add_mixture(std::mem::take(&mut edge.extracted_mixture));
        //     }
        // }

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

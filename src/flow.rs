use std::collections::HashMap;

use thiserror::Error;

use crate::{
    fluid_properties::{FluidCollection, InvalidFluidId},
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

    /// Indices of edges sending fluid from this volume.
    send_edges: Vec<usize>,
    /// Indices of edges sending fluid from to volume.
    receive_edges: Vec<usize>,
}

struct FlowEdgeState {
    /// How much liquid to move, used to calculate the `limited_move_ratio`
    move_liquid_volume: f32,
    /// How much gas to move, used to calculate the `limited_move_ratio`
    move_gas_work: f32,
    /// What pressure this edge would leave its "from" side at.
    from_pressure: f32,
    /// What pressure this edge would leave its "to" side at.
    to_pressure: f32,
    /// What proportion of the source volume is considered 100% after limiting.
    move_ratio: f32,
    /// What ratio of this edge will get moved after it has been limited by its send and receive side.
    limited_move_ratio: f32,
    extracted_mixture: Mixture,
}

/// Trait methods that allow a [FlowSimulationState] to query and modify [FluidVolume] states.
///
/// - `K` is the key that indentifies a volume.
/// - `E` is an error type that gets bubbled up through [FlowSimulationState::step] if any of the methods fail.
pub trait SimulationQueries<K, E> {
    /// Should return an iterator of all other edges that fluid can flow into from a volume.
    ///
    /// Should always return the exact same edges in the same order for any volume within one simulation step.
    fn get_edges(&self, volume_id: &K) -> Result<impl IntoIterator<Item = K>, E>;

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
        self.edges.clear();

        // Get volumes, calculate derived properties and clear working values.
        for volume_id in &self.volumes {
            let volume = queries
                .get_volume(volume_id)
                .map_err(|err| FlowQueryError(err))?;

            let state = self.volume_states.get_mut(volume_id).unwrap();

            state.volume = volume.volume();
            state.computed_properties = volume.calculate_properties(collection)?;

            state.send_edges.clear();
            state.receive_edges.clear();
        }

        // Pass over all edges and determine what amount of fluid each of them would move in isolation to equalize its two connected volumes.
        for volume_id in &self.volumes {
            let &FlowVolumeNodeSimulationState {
                volume,
                computed_properties,
                ..
            } = self.volume_states.get(volume_id).unwrap();

            if computed_properties.pressure <= 0. {
                // Zero gas pressure means no flow, skip to avoid divisions by zero.
                continue;
            }

            for other_volume_id in queries
                .get_edges(volume_id)
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

                let gas_volume = volume - computed_properties.mixture_properties.liquid_volume;
                let other_gas_volume =
                    other_volume - other_computed_properties.mixture_properties.liquid_volume;

                // What ratio of fluid from this volume needs to be moved into the other volume to equalize pressures.
                let move_ratio = {
                    let vs = volume;
                    let lt = other_computed_properties.mixture_properties.liquid_volume;
                    let ps = computed_properties.pressure;
                    let pt = other_computed_properties.pressure;

                    ((ps - pt) * gas_volume * other_gas_volume)
                        / ((ps * gas_volume * (other_gas_volume + vs))
                            + (pt * lt * other_gas_volume))
                };

                if move_ratio <= 0. || !move_ratio.is_finite() {
                    continue;
                }

                let move_liquid_volume =
                    computed_properties.mixture_properties.liquid_volume * move_ratio;
                let move_gas_work = computed_properties.pressure * gas_volume * move_ratio;

                let equalized_pressure = (computed_properties.pressure * gas_volume
                    - move_gas_work)
                    / (gas_volume + move_liquid_volume);

                println!(
                    "Moving {} to get pressure {}",
                    move_ratio, equalized_pressure
                );
                println!(
                    "will move liquid {} and work {}",
                    move_liquid_volume, move_gas_work
                );

                let edge_index = self.edges.len();

                self.edges.push(FlowEdgeState {
                    move_liquid_volume,
                    move_gas_work,
                    from_pressure: equalized_pressure,
                    to_pressure: equalized_pressure,
                    move_ratio,
                    limited_move_ratio: 0., // Will be set later.
                    extracted_mixture: Mixture::default(),
                });

                let state = self.volume_states.get_mut(volume_id).unwrap();
                state.send_edges.push(edge_index);

                let state = self.volume_states.get_mut(&other_volume_id).unwrap();
                state.receive_edges.push(edge_index);
            }
        }

        // Limit edges based on their send and receive side.

        for volume_id in &self.volumes {
            let state = self.volume_states.get(volume_id).unwrap();

            let existing_gas_work = state.computed_properties.pressure
                * (state.volume - state.computed_properties.mixture_properties.liquid_volume);

            // Limit edges adding to this volume

            let (receive_liquid_volume, receive_gas_work, pressure_target) =
                state.receive_edges.iter().fold(
                    (0., 0., state.computed_properties.pressure),
                    |(liquid, gas, max_pressure), &edge_index| {
                        let edge = self.edges.get(edge_index).unwrap();

                        (
                            liquid + edge.move_liquid_volume,
                            gas + edge.move_gas_work,
                            edge.to_pressure.max(max_pressure),
                        )
                    },
                );

            // What proportion of all the edges do we use to reach to target pressure.
            let move_ratio = (pressure_target
                * (state.volume - state.computed_properties.mixture_properties.liquid_volume)
                - existing_gas_work)
                / (receive_gas_work + pressure_target * receive_liquid_volume);

            for &edge_index in &state.receive_edges {
                let edge = self.edges.get_mut(edge_index).unwrap();

                if move_ratio.is_finite() {
                    edge.limited_move_ratio = move_ratio;

                    println!("Limiting ratio to {} from recv", move_ratio);
                } else {
                    edge.limited_move_ratio = 0.;
                }
            }

            // Limit edges taking from this volume

            let (send_liquid_volume, send_gas_work, pressure_target) =
                state.send_edges.iter().fold(
                    (0., 0., state.computed_properties.pressure),
                    |(liquid, gas, min_pressure), &edge_index| {
                        let edge = self.edges.get(edge_index).unwrap();

                        (
                            liquid + edge.move_liquid_volume,
                            gas + edge.move_gas_work,
                            edge.from_pressure.min(min_pressure),
                        )
                    },
                );

            // What proportion of all the edges do we use to reach to target pressure.
            let move_ratio = -(pressure_target
                * (state.volume - state.computed_properties.mixture_properties.liquid_volume)
                - existing_gas_work)
                / (send_gas_work + pressure_target * send_liquid_volume);

            for &edge_index in &state.send_edges {
                let edge = self.edges.get_mut(edge_index).unwrap();

                if move_ratio.is_finite() {
                    edge.limited_move_ratio = edge.limited_move_ratio.min(move_ratio);

                    println!("Limiting ratio to {} from send", move_ratio);
                } else {
                    edge.limited_move_ratio = 0.;
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

                println!("Extracting {}", edge.move_ratio * edge.limited_move_ratio);

                edge.extracted_mixture =
                    volume.mixture.extract_fluids(collection, |fluid, _| {
                        fluid.moles * edge.move_ratio * edge.limited_move_ratio
                    })?;

                dbg!(
                    volume
                        .calculate_properties(collection)
                        .unwrap()
                        .mixture_properties
                        .temperature
                );
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

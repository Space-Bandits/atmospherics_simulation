use std::collections::HashMap;

use thiserror::Error;

use crate::{
    flow::FlowSimulationError::InvalidVolumeId,
    fluid_properties::{FluidCollection, InvalidFluidId},
    fluid_volume::{ComputedVolumeProperties, FluidVolume},
    mixture::Mixture,
};

pub struct FlowSimulationState<K> {
    volumes: Vec<K>,
    volume_states: HashMap<K, FlowVolumeNodeSimulationState>,
    flow_edges: Vec<FlowEdgeState>,
    diffusion_edges: Vec<DiffusionEdgeState<K>>,
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

enum FlowEdgeState {
    Inactive,
    Active {
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
    },
}

struct DiffusionEdgeState<K> {
    volume_a_id: K,
    volume_b_id: K,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum EdgeDirectionality {
    Uni,
    Bi,
}

/// Trait methods that allow a [FlowSimulationState] to query and modify [FluidVolume] states.
///
/// - `K` is the key that indentifies a volume.
/// - `E` is an error type that gets bubbled up through [FlowSimulationState::step] if any of the methods fail.
pub trait SimulationQueries<K, E> {
    fn get_edges(&self) -> impl IntoIterator<Item = (K, K, EdgeDirectionality)>;

    fn get_volume(&self, volume_id: &K) -> Result<&FluidVolume, E>;

    fn get_volume_mut(&mut self, volume_id: &K) -> Result<&mut FluidVolume, E>;
}

impl<K> Default for FlowSimulationState<K> {
    fn default() -> Self {
        FlowSimulationState {
            volumes: Vec::new(),
            volume_states: HashMap::new(),
            flow_edges: Vec::new(),
            diffusion_edges: Vec::new(),
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
        mut queries: impl SimulationQueries<K, E>,
    ) -> Result<(), FlowSimulationError<K, E>> {
        self.flow_edges.clear();

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

        for (volume_a_id, volume_b_id, directionality) in queries.get_edges() {
            // Generate an edge always for (volume_id, other_volume_id) but also for (other_volume_id, volume_id) if bidirectional
            for (volume_id, other_volume_id) in std::iter::once((volume_a_id, volume_b_id)).chain(
                matches!(directionality, EdgeDirectionality::Bi)
                    .then_some((volume_b_id, volume_a_id)),
            ) {
                let &FlowVolumeNodeSimulationState {
                    volume,
                    computed_properties,
                    ..
                } = self
                    .volume_states
                    .get(&volume_id)
                    .ok_or(FlowSimulationError::InvalidVolumeId(volume_id))?;

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

                let edge = if move_ratio <= 0. || !move_ratio.is_finite() {
                    FlowEdgeState::Inactive
                } else {
                    let move_liquid_volume =
                        computed_properties.mixture_properties.liquid_volume * move_ratio;
                    let move_gas_work = computed_properties.pressure * gas_volume * move_ratio;

                    let equalized_pressure = (computed_properties.pressure * gas_volume
                        - move_gas_work)
                        / (gas_volume + move_liquid_volume);

                    FlowEdgeState::Active {
                        move_liquid_volume,
                        move_gas_work,
                        from_pressure: equalized_pressure,
                        to_pressure: equalized_pressure,
                        move_ratio,
                        limited_move_ratio: 0., // Will be set later.
                        extracted_mixture: Mixture::default(),
                    }
                };

                let edge_index = self.flow_edges.len();
                self.flow_edges.push(edge);

                let state = self.volume_states.get_mut(&volume_id).unwrap();
                state.send_edges.push(edge_index);

                let state = self.volume_states.get_mut(&other_volume_id).unwrap();
                state.receive_edges.push(edge_index);
            }

            if let EdgeDirectionality::Bi = directionality {
                self.diffusion_edges.push(DiffusionEdgeState {
                    volume_a_id,
                    volume_b_id,
                });
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
                        let &FlowEdgeState::Active {
                            move_liquid_volume,
                            move_gas_work,
                            to_pressure,
                            ..
                        } = self.flow_edges.get(edge_index).unwrap()
                        else {
                            return (liquid, gas, max_pressure);
                        };

                        (
                            liquid + move_liquid_volume,
                            gas + move_gas_work,
                            to_pressure.max(max_pressure),
                        )
                    },
                );

            // What proportion of all the edges do we use to reach to target pressure.
            let move_ratio = (pressure_target
                * (state.volume - state.computed_properties.mixture_properties.liquid_volume)
                - existing_gas_work)
                / (receive_gas_work + pressure_target * receive_liquid_volume);

            for &edge_index in &state.receive_edges {
                let FlowEdgeState::Active {
                    limited_move_ratio, ..
                } = self.flow_edges.get_mut(edge_index).unwrap()
                else {
                    continue;
                };

                if move_ratio.is_finite() {
                    *limited_move_ratio = move_ratio;
                } else {
                    *limited_move_ratio = 0.;
                }
            }

            // Limit edges taking from this volume

            let (send_liquid_volume, send_gas_work, pressure_target) =
                state.send_edges.iter().fold(
                    (0., 0., state.computed_properties.pressure),
                    |(liquid, gas, min_pressure), &edge_index| {
                        let &FlowEdgeState::Active {
                            move_liquid_volume,
                            move_gas_work,
                            from_pressure,
                            ..
                        } = self.flow_edges.get(edge_index).unwrap()
                        else {
                            return (liquid, gas, min_pressure);
                        };

                        (
                            liquid + move_liquid_volume,
                            gas + move_gas_work,
                            from_pressure.min(min_pressure),
                        )
                    },
                );

            // What proportion of all the edges do we use to reach to target pressure.
            let move_ratio = -(pressure_target
                * (state.volume - state.computed_properties.mixture_properties.liquid_volume)
                - existing_gas_work)
                / (send_gas_work + pressure_target * send_liquid_volume);

            for &edge_index in &state.send_edges {
                let FlowEdgeState::Active {
                    limited_move_ratio, ..
                } = self.flow_edges.get_mut(edge_index).unwrap()
                else {
                    continue;
                };

                if move_ratio.is_finite() {
                    *limited_move_ratio = limited_move_ratio.min(move_ratio);
                } else {
                    *limited_move_ratio = 0.;
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
                let FlowEdgeState::Active {
                    ref mut extracted_mixture,
                    move_ratio,
                    limited_move_ratio,
                    ..
                } = *self.flow_edges.get_mut(edge_index).unwrap()
                else {
                    continue;
                };

                *extracted_mixture = volume.mixture.extract_fluids(collection, |fluid, _| {
                    fluid.moles * move_ratio * limited_move_ratio
                })?;
            }
        }

        // Insert mixtures

        for volume_id in &self.volumes {
            let state = self.volume_states.get_mut(volume_id).unwrap();

            let volume = queries
                .get_volume_mut(volume_id)
                .map_err(|err| FlowQueryError(err))?;

            for &edge_index in &state.receive_edges {
                let FlowEdgeState::Active {
                    extracted_mixture, ..
                } = self.flow_edges.get_mut(edge_index).unwrap()
                else {
                    continue;
                };

                volume
                    .mixture
                    .add_mixture(std::mem::take(extracted_mixture));
            }

            state.computed_properties = volume.calculate_properties(collection)?;
        }

        // Calculate diffusion

        for edge in &self.diffusion_edges {
            let volume_a = queries
                .get_volume_mut(&edge.volume_a_id)
                .map_err(|err| FlowQueryError(err))?;

            let volume_b = queries
                .get_volume_mut(&edge.volume_b_id)
                .map_err(|err| FlowQueryError(err))?;

            let volume_a_properties = self
                .volume_states
                .get(&edge.volume_a_id)
                .ok_or(FlowSimulationError::InvalidVolumeId(edge.volume_a_id))?;

            let volume_b_properties = self
                .volume_states
                .get(&edge.volume_b_id)
                .ok_or(FlowSimulationError::InvalidVolumeId(edge.volume_b_id))?;
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

/// Checks if `a` and `b` have common elements assuming they are sorted.
fn does_common_edge_exist(mut a: &[usize], mut b: &[usize]) -> bool {
    loop {
        let Some(a_last) = a.last() else {
            break false;
        };

        let Some(b_last) = b.last() else {
            break false;
        };

        match a_last.cmp(b_last) {
            std::cmp::Ordering::Equal => break true,
            // if `a_last` is smaller than `b_last` then shrink `b` to get a smaller element
            std::cmp::Ordering::Less => b = &b[0..b.len() - 1],
            // if `a_last` is greater than `b_last` then shrink `a` to get a smaller element
            std::cmp::Ordering::Greater => a = &a[0..a.len() - 1],
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::flow::does_common_edge_exist;

    #[test]
    fn test_does_common_edge_exist() {
        assert!(does_common_edge_exist(&[1, 2, 3], &[1, 2, 3]));
        assert!(does_common_edge_exist(&[1, 2, 3], &[2, 3]));
        assert!(does_common_edge_exist(&[1, 2, 3], &[1, 3]));
        assert!(does_common_edge_exist(&[1, 2], &[1, 2, 3]));

        assert!(!does_common_edge_exist(&[1, 2, 3], &[4, 5]));
    }
}

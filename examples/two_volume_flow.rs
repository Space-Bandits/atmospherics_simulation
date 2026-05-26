use atmospherics_simulation::{
    ToKelvin,
    flow::{EdgeFluidPriority, FlowSimulationState, SimulationQueries},
    fluid_properties::{FluidCollection, FluidProperties, FluidTypeProperties, GasFluidProperties},
    fluid_volume::FluidVolume,
    mixture::{Fluid, Mixture},
};

fn main() {
    let collection = FluidCollection::from_iter([FluidProperties {
        fluid_type: FluidTypeProperties::Gas(GasFluidProperties {}),
        heat_capactity: 1.,
    }]);

    struct SimulationData {
        volumes: Box<[FluidVolume]>,
        edges: &'static [&'static [(usize, f32, EdgeFluidPriority)]],
    }

    let mut simulation_data = SimulationData {
        volumes: [
            FluidVolume::new(
                10.,
                Mixture::from_fluid_at_temperature(
                    &collection,
                    Fluid {
                        fluid_id: 0,
                        moles: 1.,
                    },
                    20.0.to_kelvin(),
                )
                .unwrap(),
            ),
            FluidVolume::new_empty(10.),
        ]
        .into(),
        edges: &[
            &[(1, 0.5, EdgeFluidPriority::Equal)],
            &[(0, 0.5, EdgeFluidPriority::Equal)],
        ],
    };

    impl SimulationQueries<usize, ()> for &mut SimulationData {
        fn get_edges(
            &self,
            node: &usize,
        ) -> Result<impl IntoIterator<Item = (usize, f32, EdgeFluidPriority)>, ()> {
            Ok(self.edges[*node].iter().copied())
        }

        fn get_volume(&self, volume_id: &usize) -> Result<&FluidVolume, ()> {
            self.volumes.get(*volume_id).ok_or(())
        }

        fn get_volume_mut(&mut self, volume_id: &usize) -> Result<&mut FluidVolume, ()> {
            self.volumes.get_mut(*volume_id).ok_or(())
        }
    }

    let mut simulation = FlowSimulationState::default();

    for volume_id in 0..simulation_data.volumes.len() {
        simulation.add_volume(volume_id);
    }

    simulation.step(&collection, &mut simulation_data).unwrap();
}

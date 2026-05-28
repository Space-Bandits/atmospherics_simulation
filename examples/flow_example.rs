use std::sync::LazyLock;

use atmospherics_simulation::{
    ToKelvin,
    flow::{FlowSimulationState, SimulationQueries},
    fluid_properties::{
        FluidCollection, FluidProperties, FluidTypeProperties, GasFluidProperties,
        LiquidFluidProperties,
    },
    fluid_volume::FluidVolume,
    mixture::{Fluid, Mixture},
};

static COLLECTION: LazyLock<FluidCollection> = LazyLock::new(|| {
    FluidCollection::from_iter([
        FluidProperties {
            fluid_type: FluidTypeProperties::Liquid(LiquidFluidProperties { density: 1. }),
            molar_heat_capactity: 1.,
        },
        FluidProperties {
            fluid_type: FluidTypeProperties::Gas(GasFluidProperties {}),
            molar_heat_capactity: 1.,
        },
    ])
});

struct SimulationData {
    volumes: Box<[FluidVolume]>,
    edges: &'static [&'static [usize]],
}

impl SimulationQueries<usize, ()> for SimulationData {
    fn get_edges(&self, node: &usize) -> Result<impl IntoIterator<Item = usize>, ()> {
        Ok(self.edges[*node].iter().copied())
    }

    fn get_volume(&self, volume_id: &usize) -> Result<&FluidVolume, ()> {
        self.volumes.get(*volume_id).ok_or(())
    }

    fn get_volume_mut(&mut self, volume_id: &usize) -> Result<&mut FluidVolume, ()> {
        self.volumes.get_mut(*volume_id).ok_or(())
    }
}

impl std::fmt::Debug for SimulationData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Volumes:")?;

        for volume in &self.volumes {
            let properties = volume.calculate_properties(&COLLECTION).unwrap();

            writeln!(
                f,
                "Liquid: {}, Pressure: {}",
                properties.mixture_properties.liquid_volume, properties.pressure
            )?;
        }

        Ok(())
    }
}

fn main() {
    let mut simulation_data = SimulationData {
        volumes: [
            FluidVolume::new(
                10.,
                Mixture::from_fluids_at_temperature(
                    &COLLECTION,
                    [
                        Fluid {
                            fluid_id: 0,
                            moles: 1.,
                        },
                        Fluid {
                            fluid_id: 1,
                            moles: 10.,
                        },
                    ],
                    20.0.to_kelvin(),
                )
                .unwrap(),
            ),
            FluidVolume::new(
                10.,
                Mixture::from_fluid_at_temperature(
                    &COLLECTION,
                    Fluid {
                        fluid_id: 0,
                        moles: 1.,
                    },
                    20.0.to_kelvin(),
                )
                .unwrap(),
            ),
        ]
        .into(),
        edges: &[&[1], &[0]],
    };

    let mut simulation = FlowSimulationState::default();

    for volume_id in 0..simulation_data.volumes.len() {
        simulation.add_volume(volume_id);
    }

    println!("{:?}", simulation_data);

    for _ in 0..5 {
        simulation.step(&COLLECTION, &mut simulation_data).unwrap();

        println!("{:?}", simulation_data);
    }
}

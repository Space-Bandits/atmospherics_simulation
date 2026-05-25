/// An amount of a specific fluid
#[derive(Clone, Copy, Debug)]
pub struct Fluid {
    pub fluid_id: usize,
    pub amount: f32,
}

#[derive(Default)]
pub struct Mixture {
    fluids: Vec<Fluid>,
    energy: f32,
}

impl Mixture {
    pub fn add_fluid(&mut self, fluid: Fluid) {
        match self
            .fluids
            .binary_search_by(|probe| probe.fluid_id.cmp(&fluid.fluid_id))
        {
            Ok(index) => self.fluids[index].amount += fluid.amount,
            Err(index) => self.fluids.insert(index, fluid),
        }
    }

    pub fn add_mixture(&mut self, mixture: Self) {
        for fluid in mixture.fluids {
            self.add_fluid(fluid);
        }
    }
}

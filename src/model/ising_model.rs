use crate::{graph::Graph, model::StatisticalModel};
use rand::RngExt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryCondition {
    Free,
    FixedPlus,
    FixedMinus,
}

/// Simple ferromagnetic Ising model with uniform coupling and external field.
pub struct IsingModel {
    pub coupling: f64,
    pub field: f64,
    pub temperature: f64,
    pub boundary: BoundaryCondition,
}

impl IsingModel {
    pub fn new(coupling: f64, field: f64, temperature: f64, boundary: BoundaryCondition) -> Self {
        assert!(temperature > 0.0, "Temperature must be positive");
        Self {
            coupling,
            field,
            temperature,
            boundary,
        }
    }

    fn spin_value(&self, spin: bool) -> f64 {
        if spin { 1.0 } else { -1.0 }
    }

    fn boundary_spin(&self, graph: &crate::graph::SquareLattice, site: usize) -> Option<bool> {
        if !graph.is_boundary(site) {
            return None;
        }

        match self.boundary {
            BoundaryCondition::Free => None,
            BoundaryCondition::FixedPlus => Some(true),
            BoundaryCondition::FixedMinus => Some(false),
        }
    }
}
impl StatisticalModel<crate::graph::SquareLattice> for IsingModel {
    type State = bool;

    fn initial_configuration<R: RngExt + ?Sized>(
        &self,
        graph: &crate::graph::SquareLattice,
        rng: &mut R,
    ) -> Vec<Self::State> {
        (0..graph.site_count())
            .map(|site| {
                self.boundary_spin(graph, site)
                    .unwrap_or_else(|| rng.random_bool(0.5))
            })
            .collect()
    }

    fn is_site_fixed(&self, graph: &crate::graph::SquareLattice, site: usize) -> bool {
        self.boundary_spin(graph, site).is_some()
    }

    fn energy(&self, graph: &crate::graph::SquareLattice, configuration: &[Self::State]) -> f64 {
        let mut energy = 0.0;

        for site in 0..graph.site_count() {
            let spin = self.spin_value(configuration[site]);
            energy += -self.field * spin;

            for neighbor in graph.neighbors(site) {
                if neighbor > site {
                    let neighbor_spin = self.spin_value(configuration[neighbor]);
                    energy += -self.coupling * spin * neighbor_spin;
                }
            }
        }

        energy
    }

    fn delta_energy(
        &self,
        graph: &crate::graph::SquareLattice,
        configuration: &[Self::State],
        site: usize,
        candidate: &Self::State,
    ) -> f64 {
        let current_spin = self.spin_value(configuration[site]);
        let candidate_spin = self.spin_value(*candidate);
        if current_spin == candidate_spin {
            return 0.0;
        }

        let mut delta = 0.0;
        for neighbor in graph.neighbors(site) {
            let neighbor_spin = self.spin_value(configuration[neighbor]);
            delta += -self.coupling * (candidate_spin - current_spin) * neighbor_spin;
        }
        delta += -self.field * (candidate_spin - current_spin);
        delta
    }

    fn propose_state<R: RngExt + ?Sized>(
        &self,
        current: &Self::State,
        _rng: &mut R,
    ) -> Self::State {
        !*current
    }

    fn temperature(&self) -> f64 {
        self.temperature
    }
}

#[cfg(test)]
mod tests {
    use super::{BoundaryCondition, IsingModel};
    use crate::{
        graph::{Graph, SquareLattice},
        model::StatisticalModel,
    };
    use rand::{SeedableRng, rngs::StdRng};

    #[test]
    fn fixed_plus_boundary_is_initialized_on_every_edge_site() {
        let graph = SquareLattice::new(4, 3, false);
        let model = IsingModel::new(1.0, 0.0, 2.0, BoundaryCondition::FixedPlus);
        let mut rng = StdRng::seed_from_u64(1);
        let configuration = model.initial_configuration(&graph, &mut rng);

        for site in 0..graph.site_count() {
            if graph.is_boundary(site) {
                assert!(configuration[site]);
            }
        }
    }

    #[test]
    fn fixed_minus_boundary_is_not_updated_by_sampler_hook() {
        let graph = SquareLattice::new(4, 3, false);
        let model = IsingModel::new(1.0, 0.0, 2.0, BoundaryCondition::FixedMinus);

        for site in 0..graph.site_count() {
            assert_eq!(model.is_site_fixed(&graph, site), graph.is_boundary(site));
        }
    }
}

use crate::{
    graph::{Graph, SquareLattice},
    model::StatisticalModel,
};
use rand::{Rng, RngExt};

/// q-state Potts model on a square lattice.
#[derive(Debug, Clone)]
pub struct PottsModel {
    pub q: usize,
    pub coupling: f64,
    pub temperature: f64,
}

impl PottsModel {
    pub fn new(q: usize, coupling: f64, temperature: f64) -> Self {
        assert!(q >= 2, "Potts model requires at least 2 states");
        assert!(temperature > 0.0, "Temperature must be positive");
        Self {
            q,
            coupling,
            temperature,
        }
    }
}

impl StatisticalModel<SquareLattice> for PottsModel {
    type State = u8;

    fn initial_configuration<R: Rng + ?Sized>(
        &self,
        graph: &SquareLattice,
        rng: &mut R,
    ) -> Vec<Self::State> {
        (0..graph.site_count())
            .map(|_| rng.random_range(0..self.q as u8))
            .collect()
    }

    fn energy(&self, graph: &SquareLattice, configuration: &[Self::State]) -> f64 {
        let mut energy = 0.0;
        for site in 0..graph.site_count() {
            let state = configuration[site];
            for neighbor in graph.neighbors(site) {
                if neighbor > site && configuration[neighbor] == state {
                    energy -= self.coupling;
                }
            }
        }
        energy
    }

    fn delta_energy(
        &self,
        graph: &SquareLattice,
        configuration: &[Self::State],
        site: usize,
        candidate: &Self::State,
    ) -> f64 {
        let current = configuration[site];
        if current == *candidate {
            return 0.0;
        }

        let mut old_matches = 0usize;
        let mut new_matches = 0usize;
        for neighbor in graph.neighbors(site) {
            if configuration[neighbor] == current {
                old_matches += 1;
            }
            if configuration[neighbor] == *candidate {
                new_matches += 1;
            }
        }

        self.coupling * (new_matches as f64 - old_matches as f64)
    }

    fn propose_state<R: RngExt + ?Sized>(
        &self,
        _current: &Self::State,
        rng: &mut R,
    ) -> Self::State {
        rng.random_range(0..self.q as u8)
    }

    fn temperature(&self) -> f64 {
        self.temperature
    }
}

#[cfg(test)]
mod tests {
    use super::PottsModel;
    use crate::graph::SquareLattice;
    use crate::model::StatisticalModel;

    #[test]
    fn energy_counts_matching_neighbors() {
        let lattice = SquareLattice::new(2, 2, false);
        let model = PottsModel::new(3, 1.0, 1.0);
        let config = vec![0, 0, 0, 0];
        let energy = model.energy(&lattice, &config);
        assert_eq!(energy, -4.0);
    }
}

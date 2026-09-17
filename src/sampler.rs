use crate::{
    graph::{Graph, SquareLattice},
    model::{StatisticalModel, ising_model::IsingModel, potts_model::PottsModel},
};
use rand::{RngExt, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

pub struct MetropolisSampler<'a, G: Graph, M: StatisticalModel<G>> {
    graph: &'a G,
    model: &'a M,
    configuration: Vec<M::State>,
    rng: StdRng,
}

impl<'a, G: Graph, M: StatisticalModel<G>> MetropolisSampler<'a, G, M> {
    pub fn new(graph: &'a G, model: &'a M, seed: u64) -> Self {
        let mut rng = StdRng::seed_from_u64(seed);
        let configuration = model.initial_configuration(graph, &mut rng);
        Self {
            graph,
            model,
            configuration,
            rng,
        }
    }

    pub fn run(&mut self, sweeps: usize) {
        for _ in 0..sweeps {
            for site in 0..self.graph.site_count() {
                if self.model.is_site_fixed(self.graph, site) {
                    continue;
                }
                let candidate = self
                    .model
                    .propose_state(&self.configuration[site], &mut self.rng);
                let delta =
                    self.model
                        .delta_energy(self.graph, &self.configuration, site, &candidate);
                let prob = self.model.acceptance_probability(delta);
                if self.rng.random_bool(prob) {
                    self.configuration[site] = candidate;
                }
            }
        }
    }

    pub fn configuration(&self) -> &[M::State] {
        &self.configuration
    }
}

pub struct ParallelIsingSampler<'a> {
    graph: &'a SquareLattice,
    model: &'a IsingModel,
    configuration: Vec<bool>,
    seed: u64,
}

impl<'a> ParallelIsingSampler<'a> {
    pub fn new(graph: &'a SquareLattice, model: &'a IsingModel, seed: u64) -> Self {
        assert!(
            !graph.periodic() || (graph.width() % 2 == 0 && graph.height() % 2 == 0),
            "Parallel checkerboard sampling requires even dimensions with periodic boundaries"
        );
        let mut rng = StdRng::seed_from_u64(seed);
        let configuration = model.initial_configuration(graph, &mut rng);
        Self {
            graph,
            model,
            configuration,
            seed,
        }
    }

    pub fn run(&mut self, sweeps: usize) {
        let colors = [self.sites_of_color(0), self.sites_of_color(1)];

        for sweep in 0..sweeps {
            for (color_index, sites) in colors.iter().enumerate() {
                let updates: Vec<(usize, bool)> = sites
                    .par_iter()
                    .filter_map(|&site| {
                        if self.model.is_site_fixed(self.graph, site) {
                            return None;
                        }

                        let mut rng =
                            StdRng::seed_from_u64(update_seed(self.seed, sweep, color_index, site));
                        let candidate = self
                            .model
                            .propose_state(&self.configuration[site], &mut rng);
                        let delta = self.model.delta_energy(
                            self.graph,
                            &self.configuration,
                            site,
                            &candidate,
                        );
                        let probability = self.model.acceptance_probability(delta);

                        if rng.random_bool(probability) {
                            Some((site, candidate))
                        } else {
                            None
                        }
                    })
                    .collect();

                for (site, state) in updates {
                    self.configuration[site] = state;
                }
            }
        }
    }

    pub fn configuration(&self) -> &[bool] {
        &self.configuration
    }

    fn sites_of_color(&self, color: usize) -> Vec<usize> {
        (0..self.graph.site_count())
            .filter(|&site| {
                let x = site % self.graph.width();
                let y = site / self.graph.width();
                (x + y) % 2 == color
            })
            .collect()
    }
}

fn update_seed(seed: u64, sweep: usize, color: usize, site: usize) -> u64 {
    seed.wrapping_add(sweep as u64).rotate_left(17)
        ^ (color as u64).wrapping_mul(0x9E37_79B9)
        ^ (site as u64).wrapping_mul(0x85EB_CA6B)
}

#[cfg(test)]
mod tests {
    use super::ParallelIsingSampler;
    use crate::{
        graph::{Graph, SquareLattice},
        model::ising_model::{BoundaryCondition, IsingModel},
    };

    #[test]
    fn parallel_sampler_preserves_fixed_plus_boundary() {
        let graph = SquareLattice::new(8, 8, false);
        let model = IsingModel::new(1.0, 0.0, 2.0, BoundaryCondition::FixedPlus);
        let mut sampler = ParallelIsingSampler::new(&graph, &model, 42);
        sampler.run(10);

        for site in 0..graph.site_count() {
            if graph.is_boundary(site) {
                assert!(sampler.configuration()[site]);
            }
        }
    }
}

pub struct PottsHeatBathSampler<'a> {
    graph: &'a SquareLattice,
    model: &'a PottsModel,
    configuration: Vec<u8>,
    rng: StdRng,
}

impl<'a> PottsHeatBathSampler<'a> {
    pub fn new(graph: &'a SquareLattice, model: &'a PottsModel, seed: u64) -> Self {
        let mut rng = StdRng::seed_from_u64(seed);
        let configuration = model.initial_configuration(graph, &mut rng);
        Self {
            graph,
            model,
            configuration,
            rng,
        }
    }

    pub fn run(&mut self, sweeps: usize) {
        for _ in 0..sweeps {
            for site in 0..self.graph.site_count() {
                let mut weights = vec![0.0; self.model.q];
                let mut total_weight = 0.0;

                for state in 0..self.model.q {
                    let matches = self
                        .graph
                        .neighbors(site)
                        .iter()
                        .filter(|&&neighbor| self.configuration[neighbor] == state as u8)
                        .count();
                    let weight =
                        (self.model.coupling * matches as f64 / self.model.temperature).exp();
                    weights[state] = weight;
                    total_weight += weight;
                }

                let target = self.rng.random::<f64>() * total_weight;
                let mut cumulative = 0.0;
                let mut chosen = 0usize;
                for (state, weight) in weights.iter().enumerate() {
                    cumulative += *weight;
                    if cumulative >= target {
                        chosen = state;
                        break;
                    }
                }

                self.configuration[site] = chosen as u8;
            }
        }
    }

    pub fn configuration(&self) -> &[u8] {
        &self.configuration
    }
}

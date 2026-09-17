use crate::graph::Graph;
use rand::Rng;

/// A generic statistical model on a graph.
pub trait StatisticalModel<G: Graph> {
    type State: Clone;

    fn initial_configuration<R: Rng + ?Sized>(&self, graph: &G, rng: &mut R) -> Vec<Self::State>;
    fn energy(&self, graph: &G, configuration: &[Self::State]) -> f64;
    fn delta_energy(
        &self,
        graph: &G,
        configuration: &[Self::State],
        site: usize,
        candidate: &Self::State,
    ) -> f64;
    fn propose_state<R: Rng + ?Sized>(&self, current: &Self::State, rng: &mut R) -> Self::State;
    fn temperature(&self) -> f64;

    fn is_site_fixed(&self, _graph: &G, _site: usize) -> bool {
        false
    }

    fn acceptance_probability(&self, delta_energy: f64) -> f64 {
        if delta_energy <= 0.0 {
            1.0
        } else {
            (-delta_energy / self.temperature()).exp()
        }
    }
}

pub mod ising_model;
pub mod potts_model;

impl crate::graph::Graph for () {
    fn site_count(&self) -> usize {
        0
    }

    fn neighbors(&self, _site: usize) -> Vec<usize> {
        Vec::new()
    }
}

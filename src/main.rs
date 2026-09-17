use clap::{Args as ClapArgs, Parser, Subcommand, ValueEnum};
use random_lattice_models::graph::SquareLattice;
use random_lattice_models::model::{
    StatisticalModel,
    ising_model::{BoundaryCondition, IsingModel},
    potts_model::PottsModel,
};
use random_lattice_models::output::{
    display_ising_configuration, display_potts_configuration, save_ising_configuration_image,
    save_potts_configuration_image,
};
use random_lattice_models::sampler::{ParallelIsingSampler, PottsHeatBathSampler};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum IsingBoundary {
    Free,
    Plus,
    Minus,
}

#[derive(Parser, Debug)]
#[command(name = "random_lattice_models")]
#[command(about = "Simulate classical lattice models on a finite square lattice")]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    Ising(IsingArgs),
    Potts(PottsArgs),
}

#[derive(ClapArgs, Debug)]
struct SimulationArgs {
    #[arg(long, default_value_t = 20usize)]
    width: usize,
    #[arg(long, default_value_t = 20usize)]
    height: usize,
    #[arg(long, short = 'T', default_value_t = 2.0f64)]
    temperature: f64,
    #[arg(long, default_value_t = 1000usize)]
    sweeps: usize,
    #[arg(long, default_value_t = 12345u64)]
    seed: u64,
    #[arg(long, help = "Output file name for the image")]
    output: Option<String>,
}

#[derive(ClapArgs, Debug)]
struct IsingArgs {
    #[command(flatten)]
    simulation: SimulationArgs,
    #[arg(long, default_value_t = 1.0)]
    coupling: f64,
    #[arg(long, default_value_t = 0.0, help = "Uniform external magnetic field")]
    field: f64,
    #[arg(long, value_enum, default_value_t = IsingBoundary::Free)]
    boundary: IsingBoundary,
}

#[derive(ClapArgs, Debug)]
struct PottsArgs {
    #[command(flatten)]
    simulation: SimulationArgs,
    #[arg(long, default_value_t = 4usize, help = "Number of states")]
    q: usize,
    #[arg(long, default_value_t = 1.0)]
    coupling: f64,
}

fn main() {
    let args = Args::parse();

    match args.command {
        Command::Ising(args) => {
            let simulation = args.simulation;
            let graph = SquareLattice::new(simulation.width, simulation.height, false);
            let boundary = match args.boundary {
                IsingBoundary::Free => BoundaryCondition::Free,
                IsingBoundary::Plus => BoundaryCondition::FixedPlus,
                IsingBoundary::Minus => BoundaryCondition::FixedMinus,
            };
            let model =
                IsingModel::new(args.coupling, args.field, simulation.temperature, boundary);
            let mut sampler = ParallelIsingSampler::new(&graph, &model, simulation.seed);
            sampler.run(simulation.sweeps);

            let final_config = sampler.configuration();
            let energy = model.energy(&graph, final_config);
            println!(
                "Simulating Ising model on {}x{} lattice at temperature {} for {} sweeps",
                simulation.width, simulation.height, simulation.temperature, simulation.sweeps
            );
            println!("Final energy: {:.4}\n", energy);
            display_ising_configuration(&graph, final_config);

            let image_name = simulation.output.unwrap_or_else(|| {
                format!(
                    "ising_{}x{}_T{:.2}_h{:.2}_s{}_b{}.png",
                    simulation.width,
                    simulation.height,
                    simulation.temperature,
                    args.field,
                    simulation.sweeps,
                    match boundary {
                        BoundaryCondition::Free => "_",
                        BoundaryCondition::FixedPlus => "+",
                        BoundaryCondition::FixedMinus => "-",
                    }
                )
            });
            save_ising_configuration_image(&graph, final_config, &image_name)
                .expect("Failed to write output image");
            println!("Saved image: {}", image_name);
        }
        Command::Potts(args) => {
            let simulation = args.simulation;
            let graph = SquareLattice::new(simulation.width, simulation.height, false);
            let model = PottsModel::new(args.q, args.coupling, simulation.temperature);
            let mut sampler = PottsHeatBathSampler::new(&graph, &model, simulation.seed);
            sampler.run(simulation.sweeps);

            let final_config = sampler.configuration();
            let energy = model.energy(&graph, final_config);
            println!(
                "Simulating Potts model on {}x{} lattice at temperature {} for {} sweeps",
                simulation.width, simulation.height, simulation.temperature, simulation.sweeps
            );
            println!("Final energy: {:.4}\n", energy);
            display_potts_configuration(&graph, final_config);

            let image_name = simulation.output.unwrap_or_else(|| {
                format!(
                    "potts_{}x{}_q{}_T{:.2}_s{}.png",
                    simulation.width,
                    simulation.height,
                    args.q,
                    simulation.temperature,
                    simulation.sweeps
                )
            });
            save_potts_configuration_image(&graph, final_config, &image_name, args.q)
                .expect("Failed to write output image");
            println!("Saved image: {}", image_name);
        }
    }
}

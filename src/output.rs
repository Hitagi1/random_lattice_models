use crate::graph::SquareLattice;
use image::{GrayImage, Luma, Rgb, RgbImage, imageops::resize};
use std::path::Path;

pub fn display_ising_configuration(graph: &SquareLattice, configuration: &[bool]) {
    for y in 0..graph.height() {
        for x in 0..graph.width() {
            let index = y * graph.width() + x;
            let symbol = if configuration[index] { '#' } else { '.' };
            print!("{}", symbol);
        }
        println!();
    }
}

pub fn display_potts_configuration(graph: &SquareLattice, configuration: &[u8]) {
    let symbols = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'];
    for y in 0..graph.height() {
        for x in 0..graph.width() {
            let index = y * graph.width() + x;
            let state = configuration[index] as usize;
            let symbol = symbols.get(state).copied().unwrap_or('?');
            print!("{}", symbol);
        }
        println!();
    }
}

pub fn save_ising_configuration_image<P: AsRef<Path>>(
    graph: &SquareLattice,
    configuration: &[bool],
    path: P,
) -> image::ImageResult<()> {
    let mut image = GrayImage::new(graph.width() as u32, graph.height() as u32);

    for y in 0..graph.height() {
        for x in 0..graph.width() {
            let index = y * graph.width() + x;
            let intensity = if configuration[index] { 255 } else { 0 };
            image.put_pixel(x as u32, y as u32, Luma([intensity]));
        }
    }

    if graph.width() < 1024 && graph.height() < 1024 {
        image = resize(&image, 1024, 1024, image::imageops::FilterType::Nearest);
    }

    image.save(path)
}

pub fn save_potts_configuration_image<P: AsRef<Path>>(
    graph: &SquareLattice,
    configuration: &[u8],
    path: P,
    q: usize,
) -> image::ImageResult<()> {
    let mut image = RgbImage::new(graph.width() as u32, graph.height() as u32);

    for y in 0..graph.height() {
        for x in 0..graph.width() {
            let index = y * graph.width() + x;
            let state = configuration[index] as usize;
            let color = potts_color(state, q);
            image.put_pixel(x as u32, y as u32, Rgb(color));
        }
    }

    if graph.width() < 1024 && graph.height() < 1024 {
        image = resize(&image, 1024, 1024, image::imageops::FilterType::Nearest);
    }

    image.save(path)
}

fn potts_color(state: usize, q: usize) -> [u8; 3] {
    let hue = (state as f32 * 360.0 / q.max(1) as f32) % 360.0;
    let c = 0.8;
    let x = c * (1.0 - ((hue / 60.0) % 2.0 - 1.0).abs());
    let m = 0.2;

    let (r, g, b) = match hue {
        h if h < 60.0 => (c, x, 0.0),
        h if h < 120.0 => (x, c, 0.0),
        h if h < 180.0 => (0.0, c, x),
        h if h < 240.0 => (0.0, x, c),
        h if h < 300.0 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };

    [
        ((r + m) * 255.0) as u8,
        ((g + m) * 255.0) as u8,
        ((b + m) * 255.0) as u8,
    ]
}

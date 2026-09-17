pub mod graph {
    /// A small generic graph abstraction for finite site-based models.
    pub trait Graph {
        fn site_count(&self) -> usize;
        fn neighbors(&self, site: usize) -> Vec<usize>;
    }

    /// A finite subsquare of Z^2 with optional periodic boundary conditions.
    #[derive(Debug, Clone, Copy)]
    pub struct SquareLattice {
        width: usize,
        height: usize,
        periodic: bool,
    }

    impl SquareLattice {
        pub fn new(width: usize, height: usize, periodic: bool) -> Self {
            assert!(
                width > 0 && height > 0,
                "Lattice dimensions must be positive"
            );
            Self {
                width,
                height,
                periodic,
            }
        }

        pub fn width(&self) -> usize {
            self.width
        }

        pub fn height(&self) -> usize {
            self.height
        }

        pub fn periodic(&self) -> bool {
            self.periodic
        }

        pub fn is_boundary(&self, site: usize) -> bool {
            let x = site % self.width;
            let y = site / self.width;
            !self.periodic && (x == 0 || x + 1 == self.width || y == 0 || y + 1 == self.height)
        }

        fn index(&self, x: usize, y: usize) -> usize {
            y * self.width + x
        }

        fn in_bounds(&self, x: isize, width: usize) -> bool {
            x >= 0 && (x as usize) < width
        }
    }

    impl Graph for SquareLattice {
        fn site_count(&self) -> usize {
            self.width * self.height
        }

        fn neighbors(&self, site: usize) -> Vec<usize> {
            let x = site % self.width;
            let y = site / self.width;
            let mut result = Vec::with_capacity(4);

            let candidates = [
                (x as isize - 1, y as isize),
                (x as isize + 1, y as isize),
                (x as isize, y as isize - 1),
                (x as isize, y as isize + 1),
            ];

            for (nx, ny) in candidates {
                if self.periodic {
                    let wrapped_x =
                        ((nx % self.width as isize) + self.width as isize) % self.width as isize;
                    let wrapped_y =
                        ((ny % self.height as isize) + self.height as isize) % self.height as isize;
                    result.push(self.index(wrapped_x as usize, wrapped_y as usize));
                } else if self.in_bounds(nx, self.width) && self.in_bounds(ny, self.height) {
                    result.push(self.index(nx as usize, ny as usize));
                }
            }

            result
        }
    }
}

pub mod model;

pub mod sampler;

pub mod output;

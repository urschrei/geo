use robust::Coord;

use crate::{Coord as GeoCoord, GeoNum};

/// Utility functions for working with quadrants of the cartesian plane,
/// which are labeled as follows:
/// ```ignore
///          (+)
///        NW ┃ NE
///    (-) ━━━╋━━━━ (+)
///        SW ┃ SE
///          (-)
/// ```
#[derive(Debug, Clone, Copy, PartialOrd, PartialEq, Eq)]
pub enum Quadrant {
    NE,
    NW,
    SW,
    SE,
}

impl Quadrant {
    pub fn new<F: GeoNum>(dx: F, dy: F) -> Option<Quadrant> {
        if dx.is_zero() && dy.is_zero() {
            return None;
        }

        match (dy >= F::zero(), dx >= F::zero()) {
            (true, true) => Quadrant::NE,
            (true, false) => Quadrant::NW,
            (false, false) => Quadrant::SW,
            (false, true) => Quadrant::SE,
        }
        .into()
    }
    pub fn new_from_coords<F: GeoNum>(c0: GeoCoord<F>, c1: GeoCoord<F>) -> Option<Quadrant> {
        if c1.x == c0.x && c1.y == c0.y {
            return None;
        }
        if c1.x >= c0.x {
            if c1.y >= c0.y {
                Some(Quadrant::NE)
            } else {
                Some(Quadrant::SE)
            }
        } else if c1.y >= c0.y {
            Some(Quadrant::NW)
        } else {
            Some(Quadrant::SW)
        }
    }
}

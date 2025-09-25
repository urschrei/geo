use super::{Distance, Euclidean};
use crate::HasDimensions;
use crate::algorithm::BoundingRect;
use crate::algorithm::Intersects;
use crate::algorithm::centroid::Centroid;
use crate::algorithm::coords_iter::CoordsIter;
use crate::algorithm::vector_ops::Vector2DOps;
use crate::coordinate_position::{CoordPos, coord_pos_relative_to_ring};
use crate::geometry::*;
use crate::{CoordFloat, GeoFloat, GeoNum};
use num_traits::{Bounded, Float};
use rstar::RTree;
use rstar::primitives::CachedEnvelope;

// Distance is a symmetric operation, so we can implement it once for both
macro_rules! symmetric_distance_impl {
    ($t:ident, $a:ty, $b:ty) => {
        impl<F> $crate::Distance<F, $a, $b> for Euclidean
        where
            F: $t,
        {
            fn distance(&self, a: $a, b: $b) -> F {
                self.distance(b, a)
            }
        }
    };
}

// ┌───────────────────────────┐
// │ Implementations for Coord │
// └───────────────────────────┘

impl<F: CoordFloat> Distance<F, Coord<F>, Coord<F>> for Euclidean {
    fn distance(&self, origin: Coord<F>, destination: Coord<F>) -> F {
        let delta = origin - destination;
        delta.x.hypot(delta.y)
    }
}
impl<F: CoordFloat> Distance<F, Coord<F>, &Line<F>> for Euclidean {
    fn distance(&self, coord: Coord<F>, line: &Line<F>) -> F {
        ::geo_types::private_utils::point_line_euclidean_distance(Point(coord), *line)
    }
}

// ┌───────────────────────────┐
// │ Implementations for Point │
// └───────────────────────────┘

/// Calculate the Euclidean distance (a.k.a. pythagorean distance) between two Points
impl<F: CoordFloat> Distance<F, Point<F>, Point<F>> for Euclidean {
    /// Calculate the Euclidean distance (a.k.a. pythagorean distance) between two Points
    ///
    /// # Units
    /// - `origin`, `destination`: Point where the units of x/y represent non-angular units,
    ///   e.g. meters or miles, not lon/lat. For lon/lat points, use the
    ///   [`Haversine`] or [`Geodesic`] [metric spaces].
    /// - returns: distance in the same units as the `origin` and `destination` points
    ///
    /// # Example
    /// ```
    /// use geo::{Euclidean, Distance};
    /// use geo::Point;
    /// // web mercator
    /// let new_york_city = Point::new(-8238310.24, 4942194.78);
    /// // web mercator
    /// let london = Point::new(-14226.63, 6678077.70);
    /// let distance: f64 = Euclidean.distance(new_york_city, london);
    ///
    /// assert_eq!(
    ///     8_405_286., // meters in web mercator
    ///     distance.round()
    /// );
    /// ```
    ///
    /// [`Haversine`]: crate::line_measures::metric_spaces::Haversine
    /// [`Geodesic`]: crate::line_measures::metric_spaces::Geodesic
    /// [metric spaces]: crate::line_measures::metric_spaces
    fn distance(&self, origin: Point<F>, destination: Point<F>) -> F {
        self.distance(origin.0, destination.0)
    }
}

impl<F: CoordFloat> Distance<F, &Point<F>, &Point<F>> for Euclidean {
    fn distance(&self, origin: &Point<F>, destination: &Point<F>) -> F {
        self.distance(*origin, *destination)
    }
}

impl<F: CoordFloat> Distance<F, &Point<F>, &Line<F>> for Euclidean {
    fn distance(&self, origin: &Point<F>, destination: &Line<F>) -> F {
        geo_types::private_utils::point_line_euclidean_distance(*origin, *destination)
    }
}

impl<F: CoordFloat> Distance<F, &Point<F>, &LineString<F>> for Euclidean {
    fn distance(&self, origin: &Point<F>, destination: &LineString<F>) -> F {
        geo_types::private_utils::point_line_string_euclidean_distance(*origin, destination)
    }
}

impl<F: GeoFloat> Distance<F, &Point<F>, &Polygon<F>> for Euclidean {
    fn distance(&self, point: &Point<F>, polygon: &Polygon<F>) -> F {
        // No need to continue if the polygon intersects the point, or is zero-length
        if polygon.exterior().0.is_empty() || polygon.intersects(point) {
            return F::zero();
        }
        // fold the minimum interior ring distance if any, followed by the exterior
        // shell distance, returning the minimum of the two distances
        polygon
            .interiors()
            .iter()
            .map(|ring| self.distance(point, ring))
            .fold(Bounded::max_value(), |accum: F, val| accum.min(val))
            .min(
                polygon
                    .exterior()
                    .lines()
                    .map(|line| {
                        ::geo_types::private_utils::line_segment_distance(
                            point.0, line.start, line.end,
                        )
                    })
                    .fold(Bounded::max_value(), |accum, val| accum.min(val)),
            )
    }
}

// ┌──────────────────────────┐
// │ Implementations for Line │
// └──────────────────────────┘

symmetric_distance_impl!(CoordFloat, &Line<F>, Coord<F>);
symmetric_distance_impl!(CoordFloat, &Line<F>, &Point<F>);

impl<F: GeoFloat> Distance<F, &Line<F>, &Line<F>> for Euclidean {
    fn distance(&self, line_a: &Line<F>, line_b: &Line<F>) -> F {
        if line_a.intersects(line_b) {
            return F::zero();
        }
        // minimum of all Point-Line distances
        self.distance(&line_a.start_point(), line_b)
            .min(self.distance(&line_a.end_point(), line_b))
            .min(self.distance(&line_b.start_point(), line_a))
            .min(self.distance(&line_b.end_point(), line_a))
    }
}

impl<F: GeoFloat> Distance<F, &Line<F>, &LineString<F>> for Euclidean {
    fn distance(&self, line: &Line<F>, line_string: &LineString<F>) -> F {
        line_string
            .lines()
            .fold(Bounded::max_value(), |acc, segment| {
                acc.min(self.distance(line, &segment))
            })
    }
}

impl<F: GeoFloat> Distance<F, &Line<F>, &Polygon<F>> for Euclidean {
    fn distance(&self, line: &Line<F>, polygon: &Polygon<F>) -> F {
        if line.intersects(polygon) {
            return F::zero();
        }

        // REVIEW: This impl changed slightly.
        std::iter::once(polygon.exterior())
            .chain(polygon.interiors().iter())
            .fold(Bounded::max_value(), |acc, line_string| {
                acc.min(self.distance(line, line_string))
            })
    }
}

// ┌────────────────────────────────┐
// │ Implementations for LineString │
// └────────────────────────────────┘

symmetric_distance_impl!(CoordFloat, &LineString<F>, &Point<F>);
symmetric_distance_impl!(GeoFloat, &LineString<F>, &Line<F>);

impl<F: GeoFloat> Distance<F, &LineString<F>, &LineString<F>> for Euclidean {
    fn distance(&self, line_string_a: &LineString<F>, line_string_b: &LineString<F>) -> F {
        if line_string_a.intersects(line_string_b) {
            F::zero()
        } else {
            nearest_neighbour_distance(line_string_a, line_string_b)
        }
    }
}

impl<F: GeoFloat> Distance<F, &LineString<F>, &Polygon<F>> for Euclidean {
    fn distance(&self, line_string: &LineString<F>, polygon: &Polygon<F>) -> F {
        if line_string.intersects(polygon) {
            F::zero()
        } else if !polygon.interiors().is_empty()
            // FIXME: Explodes on empty line_string
            && ring_contains_coord(polygon.exterior(), line_string.0[0])
        {
            // check each ring distance, returning the minimum
            let mut mindist: F = Float::max_value();
            for ring in polygon.interiors() {
                mindist = mindist.min(nearest_neighbour_distance(line_string, ring))
            }
            mindist
        } else {
            nearest_neighbour_distance(line_string, polygon.exterior())
        }
    }
}

// ┌─────────────────────────────┐
// │ Implementations for Polygon │
// └─────────────────────────────┘

symmetric_distance_impl!(GeoFloat, &Polygon<F>, &Point<F>);
symmetric_distance_impl!(GeoFloat, &Polygon<F>, &Line<F>);
symmetric_distance_impl!(GeoFloat, &Polygon<F>, &LineString<F>);
impl<F: GeoFloat> Distance<F, &Polygon<F>, &Polygon<F>> for Euclidean {
    fn distance(&self, polygon_a: &Polygon<F>, polygon_b: &Polygon<F>) -> F {
        if polygon_a.intersects(polygon_b) {
            return F::zero();
        }
        if polygon_a.is_empty() || polygon_b.is_empty() {
            return F::zero();
        }

        // Check if bounding boxes are non-overlapping: if so, use project-and-sort optimisation
        let bbox_a = polygon_a.bounding_rect();
        let bbox_b = polygon_b.bounding_rect();

        if let (Some(rect_a), Some(rect_b)) = (bbox_a, bbox_b) {
            // Check for bbox separation along both axes
            // TODO: do we have anything built-in that does this cheaply?
            let x_separated = rect_a.max().x < rect_b.min().x || rect_b.max().x < rect_a.min().x;
            let y_separated = rect_a.max().y < rect_b.min().y || rect_b.max().y < rect_a.min().y;

            if x_separated || y_separated {
                return minimum_separable_distance(polygon_a, polygon_b);
            }
        }

        // FIXME: explodes when polygon_b.exterior() is empty
        // Containment check
        if !polygon_a.interiors().is_empty()
            && ring_contains_coord(polygon_a.exterior(), polygon_b.exterior().0[0])
        {
            // check each ring distance, returning the minimum
            let mut mindist: F = Float::max_value();
            for ring in polygon_a.interiors() {
                mindist = mindist.min(nearest_neighbour_distance(polygon_b.exterior(), ring))
            }
            return mindist;
        } else if !polygon_b.interiors().is_empty()
            // FIXME: explodes when polygon_a.exterior() is empty
            && ring_contains_coord(polygon_b.exterior(), polygon_a.exterior().0[0])
        {
            let mut mindist: F = Float::max_value();
            for ring in polygon_b.interiors() {
                mindist = mindist.min(nearest_neighbour_distance(polygon_a.exterior(), ring))
            }
            return mindist;
        }
        nearest_neighbour_distance(polygon_a.exterior(), polygon_b.exterior())
    }
}

// ┌────────────────────────────────────────┐
// │ Implementations for Rect and Triangle  │
// └────────────────────────────────────────┘

/// Implements Euclidean distance for Triangles and Rects by converting them to polygons.
macro_rules! impl_euclidean_distance_for_polygonlike_geometry {
  ($polygonlike:ty,  [$($geometry_b:ty),*]) => {
      impl<F: GeoFloat> Distance<F, $polygonlike, $polygonlike> for Euclidean
      {
          fn distance(&self, origin: $polygonlike, destination: $polygonlike) -> F {
              self.distance(&origin.to_polygon(), destination)
          }
      }
      $(
          impl<F: GeoFloat> Distance<F, $polygonlike, $geometry_b> for Euclidean
          {
              fn distance(&self, polygonlike: $polygonlike, geometry_b: $geometry_b) -> F {
                    self.distance(&polygonlike.to_polygon(), geometry_b)
              }
          }
          symmetric_distance_impl!(GeoFloat, $geometry_b, $polygonlike);
      )*
  };
}

impl_euclidean_distance_for_polygonlike_geometry!(&Triangle<F>,  [&Point<F>, &Line<F>, &LineString<F>, &Polygon<F>, &Rect<F>]);
impl_euclidean_distance_for_polygonlike_geometry!(&Rect<F>,      [&Point<F>, &Line<F>, &LineString<F>, &Polygon<F>]);

// ┌───────────────────────────────────────────┐
// │ Implementations for multi geometry types  │
// └───────────────────────────────────────────┘

/// Euclidean distance implementation for multi geometry types.
macro_rules! impl_euclidean_distance_for_iter_geometry {
    ($iter_geometry:ty,  [$($to_geometry:ty),*]) => {
        impl<F: GeoFloat> Distance<F, $iter_geometry, $iter_geometry> for Euclidean {
            fn distance(&self, origin: $iter_geometry, destination: $iter_geometry) -> F {
                origin
                    .iter()
                    .fold(Bounded::max_value(), |accum: F, member| {
                        accum.min(self.distance(member, destination))
                    })
             }
        }
        $(
            impl<F: GeoFloat> Distance<F, $iter_geometry, $to_geometry> for Euclidean {
                fn distance(&self, iter_geometry: $iter_geometry, to_geometry: $to_geometry) -> F {
                    iter_geometry
                        .iter()
                        .fold(Bounded::max_value(), |accum: F, member| {
                            accum.min(self.distance(member, to_geometry))
                        })
                }
            }
            symmetric_distance_impl!(GeoFloat, $to_geometry, $iter_geometry);
        )*
  };
}

impl_euclidean_distance_for_iter_geometry!(&MultiPoint<F>,         [&Point<F>, &Line<F>, &LineString<F>, &MultiLineString<F>, &Polygon<F>, &MultiPolygon<F>, &GeometryCollection<F>, &Rect<F>, &Triangle<F>]);
impl_euclidean_distance_for_iter_geometry!(&MultiLineString<F>,    [&Point<F>, &Line<F>, &LineString<F>,                      &Polygon<F>, &MultiPolygon<F>, &GeometryCollection<F>, &Rect<F>, &Triangle<F>]);
impl_euclidean_distance_for_iter_geometry!(&MultiPolygon<F>,       [&Point<F>, &Line<F>, &LineString<F>,                      &Polygon<F>,                   &GeometryCollection<F>, &Rect<F>, &Triangle<F>]);
impl_euclidean_distance_for_iter_geometry!(&GeometryCollection<F>, [&Point<F>, &Line<F>, &LineString<F>,                      &Polygon<F>,                                           &Rect<F>, &Triangle<F>]);

// ┌──────────────────────────────┐
// │ Implementation for Geometry  │
// └──────────────────────────────┘

/// Euclidean distance implementation for every specific Geometry type to Geometry<T>.
macro_rules! impl_euclidean_distance_for_geometry_and_variant {
  ([$($target:ty),*]) => {
      $(
          impl<F: GeoFloat> Distance<F, $target, &Geometry<F>> for Euclidean {
              fn distance(&self, origin: $target, destination: &Geometry<F>) -> F {
                  match destination {
                      Geometry::Point(point) => self.distance(origin, point),
                      Geometry::Line(line) => self.distance(origin, line),
                      Geometry::LineString(line_string) => self.distance(origin, line_string),
                      Geometry::Polygon(polygon) => self.distance(origin, polygon),
                      Geometry::MultiPoint(multi_point) => self.distance(origin, multi_point),
                      Geometry::MultiLineString(multi_line_string) => self.distance(origin, multi_line_string),
                      Geometry::MultiPolygon(multi_polygon) => self.distance(origin, multi_polygon),
                      Geometry::GeometryCollection(geometry_collection) => self.distance(origin, geometry_collection),
                      Geometry::Rect(rect) => self.distance(origin, rect),
                      Geometry::Triangle(triangle) => self.distance(origin, triangle),
                  }
              }
          }
          symmetric_distance_impl!(GeoFloat, &Geometry<F>, $target);
      )*
  };
}

impl_euclidean_distance_for_geometry_and_variant!([&Point<F>, &MultiPoint<F>, &Line<F>, &LineString<F>, &MultiLineString<F>, &Polygon<F>, &MultiPolygon<F>, &Triangle<F>, &Rect<F>, &GeometryCollection<F>]);

impl<F: GeoFloat> Distance<F, &Geometry<F>, &Geometry<F>> for Euclidean {
    fn distance(&self, origin: &Geometry<F>, destination: &Geometry<F>) -> F {
        match origin {
            Geometry::Point(point) => self.distance(point, destination),
            Geometry::Line(line) => self.distance(line, destination),
            Geometry::LineString(line_string) => self.distance(line_string, destination),
            Geometry::Polygon(polygon) => self.distance(polygon, destination),
            Geometry::MultiPoint(multi_point) => self.distance(multi_point, destination),
            Geometry::MultiLineString(multi_line_string) => {
                self.distance(multi_line_string, destination)
            }
            Geometry::MultiPolygon(multi_polygon) => self.distance(multi_polygon, destination),
            Geometry::GeometryCollection(geometry_collection) => {
                self.distance(geometry_collection, destination)
            }
            Geometry::Rect(rect) => self.distance(rect, destination),
            Geometry::Triangle(triangle) => self.distance(triangle, destination),
        }
    }
}

// ┌───────────────────────────┐
// │ Implementations utilities │
// └───────────────────────────┘

/// Uses an R* tree and nearest-neighbour lookups to calculate minimum distances
// This is somewhat slow and memory-inefficient, but certainly better than quadratic time
fn nearest_neighbour_distance<F: GeoFloat>(geom1: &LineString<F>, geom2: &LineString<F>) -> F {
    let tree_a = RTree::bulk_load(geom1.lines().map(CachedEnvelope::new).collect());
    let tree_b = RTree::bulk_load(geom2.lines().map(CachedEnvelope::new).collect());
    // Return minimum distance between all geom a points and geom b lines, and all geom b points and geom a lines
    geom2
        .points()
        .fold(Bounded::max_value(), |acc: F, point| {
            let nearest = tree_a.nearest_neighbor(&point).unwrap();
            acc.min(Euclidean.distance(nearest as &Line<F>, &point))
        })
        .min(geom1.points().fold(Bounded::max_value(), |acc, point| {
            let nearest = tree_b.nearest_neighbor(&point).unwrap();
            acc.min(Euclidean.distance(nearest as &Line<F>, &point))
        }))
}

fn ring_contains_coord<T: GeoNum>(ring: &LineString<T>, c: Coord<T>) -> bool {
    match coord_pos_relative_to_ring(c, ring) {
        CoordPos::Inside => true,
        CoordPos::OnBoundary | CoordPos::Outside => false,
    }
}

/// Projects polygon vertices onto an arbitrary axis using scalar projection.
///
/// This function implements the projection step of PostGIS's "project and sort" algorithm
/// for efficient polygon distance calculation. It uses the scalar projection formula to
/// map 2D vertex coordinates to 1D distances along the specified axis.
///
/// ## Mathematical Basis
///
/// For each vertex `v`, computes the parametric projection onto the line segment
/// defined by `axis_start` → `axis_end`. Given a line segment from point A to B
/// and a point P to project, the parametric coordinate r is:
///
/// ```text
/// r = (P - A) · (B - A) / ||B - A||²
/// ```
///
/// where:
/// - `·` denotes the dot product
/// - `||B - A||²` is the squared Euclidean distance
/// - `r = 0` corresponds to point A (axis_start)
/// - `r = 1` corresponds to point B (axis_end)
/// - `r < 0` indicates positions before A
/// - `r > 1` indicates positions beyond B
///
/// ## PostGIS Implementation
///
/// This formula is used in PostGIS's distance calculation function:
/// [`lw_dist2d_pt_seg()`](https://github.com/postgis/postgis/blob/bd63102e05c5df317cb631efe782cc79e3e053db/liblwgeom/measures.c#L2305)
/// to find the closest point on a line segment to a given point.
///
/// The parametric approach differs from orthogonal projection `(a · b) / ||b||`
/// by dividing by ||b||² instead of ||b||, yielding normalized coordinates
/// along the line segment rather than absolute distances.
///
/// ### Why Parametric Projection
///
/// The parametric value directly indicates where along the axis the projection falls,
/// making it easy to sort points by their position along the axis. This enables
/// PostGIS's "project and sort" optimization: vertices are projected onto an axis
/// perpendicular to the line between polygon centers, then sorted to efficiently
/// identify which segments are worth comparing for distance calculations.
///
/// ## Algorithm Origin
///
/// Based on PostGIS's distance optimization described in:
/// ["Inside PostGIS: Calculating Distance"](https://www.crunchydata.com/blog/inside-postgis-calculating-distance)
///
/// ## Implementation Notes
///
/// - Only projects exterior ring vertices (interior rings cannot be closest to external polygons)
/// - Uses squared axis length to avoid expensive square root calculation
/// - Returns (parametric_coordinate, original_vertex_index) for spatial indexing
fn project_vertices_onto_axis<F: GeoFloat>(
    poly: &Polygon<F>,
    axis_start: Point<F>,
    axis_end: Point<F>,
) -> Vec<(F, usize)> {
    let c1 = axis_start.0;
    let c2 = axis_end.0;

    // Axis vector from start to end
    let axis_vector = Coord {
        x: c2.x - c1.x,
        y: c2.y - c1.y,
    };

    // Pre-calculate axis length squared (avoid sqrt)
    let axis_length_squared = axis_vector.magnitude_squared();

    // Project exterior ring vertices
    poly.exterior()
        .coords_iter()
        .enumerate()
        .map(|(idx, vertex_coord)| {
            // Vector from axis start to vertex
            let vertex_vector = Coord {
                x: vertex_coord.x - c1.x,
                y: vertex_coord.y - c1.y,
            };

            let dot = vertex_vector.dot_product(axis_vector);

            // Normalize by axis length squared to get projection distance
            let projection_distance = dot / axis_length_squared;

            (projection_distance, idx)
        })
        .collect()
}

/// Calculates minimum distance between two horizontally-separated polygons using
/// PostGIS-style "project and sort" algorithm.
///
/// This approach:
/// 1. Calculates centroids of both polygons
/// 2. Projects all vertices onto the line between centroids
/// 3. Sorts projections to find candidate closest vertices
/// 4. Uses segment context to calculate accurate distance
///
/// This handles complex non-convex polygons correctly, unlike simple x-extremes approach.
fn project_and_sort_distance<F: GeoFloat>(poly1: &Polygon<F>, poly2: &Polygon<F>) -> F {
    // Step 1: Get centroids
    let c1 = poly1
        .centroid()
        .expect("Non-empty polygon should have centroid");
    let c2 = poly2
        .centroid()
        .expect("Non-empty polygon should have centroid");

    // Step 2: Project all vertices onto the centroid-to-centroid axis
    let mut proj1 = project_vertices_onto_axis(poly1, c1, c2);
    let mut proj2 = project_vertices_onto_axis(poly2, c1, c2);

    // Step 3: Sort by projection distance
    proj1.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    proj2.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // Step 4: Determine which polygon is "left" vs "right" along the projection axis
    // and get candidate vertices (rightmost of left, leftmost of right)
    let (left_poly, right_poly, left_candidate_idx, right_candidate_idx) =
        if proj1.last().unwrap().0 < proj2.first().unwrap().0 {
            // poly1 is left, poly2 is right
            (
                poly1,
                poly2,
                proj1
                    .last()
                    .expect("The 'left' Polygon projection should not be empty, but it is")
                    .1,
                proj2
                    .first()
                    .expect("The 'right' Polygon projection should not be empty, but it is")
                    .1,
            )
        } else {
            // poly2 is left, poly1 is right
            (
                poly2,
                poly1,
                proj2.last().unwrap().1,
                proj1.first().unwrap().1,
            )
        };

    // Step 5: Use segment context to ensure correct distance calculation
    let left_linestring = create_context_linestring(left_poly, left_candidate_idx);
    let right_linestring = create_context_linestring(right_poly, right_candidate_idx);

    // Step 6: Calculate final distance between the context LineStrings
    Euclidean.distance(&left_linestring, &right_linestring)
}

/// Calculates the minimum separable distance between two non-overlapping polygons.
///
/// This function uses a PostGIS-inspired "project and sort" algorithm:
/// 1. Projects both polygons onto the line between their centroids
/// 2. Sorts the projected points to find candidates
/// 3. Uses segment context around candidates for accurate distance calculation
///
/// This approach correctly handles complex non-convex polygons where simple
/// x-extremes approaches fail.
///
/// # Arguments
/// * `poly1` - First polygon (assumed non-overlapping with poly2)
/// * `poly2` - Second polygon (assumed non-overlapping with poly1)
///
/// # Returns
/// The minimum separable distance between the two polygons
pub(crate) fn minimum_separable_distance<F: GeoFloat>(poly1: &Polygon<F>, poly2: &Polygon<F>) -> F {
    project_and_sort_distance(poly1, poly2)
}

/// Creates a 3-vertex LineString from a Polygon index, with the given index as the middle vertex
fn create_context_linestring<F: GeoFloat>(poly: &Polygon<F>, center_idx: usize) -> LineString<F> {
    // Only consider the exterior ring since interior rings cannot be closest to external polygons
    let exterior_coords = &poly.exterior().0;
    let num_coords = exterior_coords.len();

    if num_coords < 3 || center_idx >= num_coords {
        // Fallback: return a degenerate LineString with just the first point if available
        if let Some(&coord) = exterior_coords.first() {
            return LineString::new(vec![coord]);
        } else {
            return LineString::new(vec![]);
        }
    }

    // For closed polygons, the last coordinate duplicates the first
    // We have n-1 unique vertices (indices 0 to n-2) plus closing duplicate at n-1
    let effective_coords =
        if num_coords > 1 && exterior_coords[0] == exterior_coords[num_coords - 1] {
            num_coords - 1 // Exclude the duplicate closing vertex
        } else {
            num_coords
        };

    // Calculate previous and next indices with wraparound for closed polygon
    let prev_idx = if center_idx == 0 {
        effective_coords.saturating_sub(1)
    } else {
        center_idx - 1
    };
    let next_idx = if center_idx >= effective_coords - 1 {
        0 // Wrap to start for closed polygon
    } else {
        center_idx + 1
    };

    // Create a 3-vertex LineString using direct indexing
    LineString::new(vec![
        exterior_coords[prev_idx],
        exterior_coords[center_idx],
        exterior_coords[next_idx],
    ])
}
#[cfg(test)]
mod test {
    use super::*;
    use crate::orient::{Direction, Orient};
    use crate::{
        Convert, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
    };
    use geo_types::wkt;
    use geo_types::{coord, polygon, private_utils::line_segment_distance};

    #[test]
    fn test_y_axis_separated_polygons() {
        // Two polygons separated along the y-axis (one above the other)
        let poly1: Polygon = wkt! { POLYGON((0 0, 10 0, 10 5, 0 5, 0 0)) }.convert();
        let poly2: Polygon = wkt! { POLYGON((2 10, 8 10, 8 15, 2 15, 2 10)) }.convert();

        // Bounding boxes overlap in x-axis but are separated in y-axis
        // This should now trigger the project-and-sort algorithm
        let distance = Euclidean.distance(&poly1, &poly2);
        assert_relative_eq!(distance, 5.0); // Distance from top of poly1 to bottom of poly2
    }

    #[test]
    fn test_separable_distance() {
        // Sanity check: `a` is a triangle pointing right. `b` is a rectangle strictly right of `a`.
        //             _
        // a->  |>    |_| <-b
        let a: Polygon = wkt! { POLYGON ((160 70, 160 320, 400 200, 160 70)) }.convert();
        let b: Polygon = wkt! { POLYGON ((500 100, 500 300, 700 300, 700 100, 500 100)) }.convert();
        assert_eq!(100.0, Euclidean.distance(&a, &b));

        // Leave the base of `b`, but extend it's upper area to trigger a pathological case in the overlap heuristic.
        // Intuitively, since the closest edge (the left edge of `b`) is still there, it can't be farther than
        // the previous test case.
        //          ____
        //         [_   | <-b
        //           \  |
        //            | |
        //            | |
        // a->  |>    |_|
        let a: Polygon = wkt! { POLYGON ((160 70, 160 320, 400 200, 160 70)) }.convert();
        let b: Polygon =
            wkt! { POLYGON ((450 500, 485 449, 500 400, 500 200, 900 500, 640 520, 450 500)) }
                .convert();
        // FAILS!
        // assertion `left == right` failed
        // left: 100.0
        // right: 260.72552617647546
        assert_eq!(100.0, Euclidean.distance(&a, &b));
    }

    #[test]
    fn line_segment_distance_test() {
        let o1 = Point::new(8.0, 0.0);
        let o2 = Point::new(5.5, 0.0);
        let o3 = Point::new(5.0, 0.0);
        let o4 = Point::new(4.5, 1.5);

        let p1 = Point::new(7.2, 2.0);
        let p2 = Point::new(6.0, 1.0);

        let dist = line_segment_distance(o1, p1, p2);
        let dist2 = line_segment_distance(o2, p1, p2);
        let dist3 = line_segment_distance(o3, p1, p2);
        let dist4 = line_segment_distance(o4, p1, p2);
        // Results agree with Shapely
        assert_relative_eq!(dist, 2.0485900789263356);
        assert_relative_eq!(dist2, 1.118033988749895);
        assert_relative_eq!(dist3, std::f64::consts::SQRT_2); // workaround clippy::correctness error approx_constant (1.4142135623730951)
        assert_relative_eq!(dist4, 1.5811388300841898);
        // Point is on the line
        let zero_dist = line_segment_distance(p1, p1, p2);
        assert_relative_eq!(zero_dist, 0.0);
    }
    #[test]
    // Point to Polygon, outside point
    fn point_polygon_distance_outside_test() {
        // an octagon
        let points = vec![
            (5., 1.),
            (4., 2.),
            (4., 3.),
            (5., 4.),
            (6., 4.),
            (7., 3.),
            (7., 2.),
            (6., 1.),
            (5., 1.),
        ];
        let ls = LineString::from(points);
        let poly = Polygon::new(ls, vec![]);
        // A Random point outside the octagon
        let p = Point::new(2.5, 0.5);
        let dist = Euclidean.distance(&p, &poly);
        assert_relative_eq!(dist, 2.1213203435596424);
    }
    #[test]
    // Point to Polygon, inside point
    fn point_polygon_distance_inside_test() {
        // an octagon
        let points = vec![
            (5., 1.),
            (4., 2.),
            (4., 3.),
            (5., 4.),
            (6., 4.),
            (7., 3.),
            (7., 2.),
            (6., 1.),
            (5., 1.),
        ];
        let ls = LineString::from(points);
        let poly = Polygon::new(ls, vec![]);
        // A Random point inside the octagon
        let p = Point::new(5.5, 2.1);
        let dist = Euclidean.distance(&p, &poly);
        assert_relative_eq!(dist, 0.0);
    }
    #[test]
    // Point to Polygon, on boundary
    fn point_polygon_distance_boundary_test() {
        // an octagon
        let points = vec![
            (5., 1.),
            (4., 2.),
            (4., 3.),
            (5., 4.),
            (6., 4.),
            (7., 3.),
            (7., 2.),
            (6., 1.),
            (5., 1.),
        ];
        let ls = LineString::from(points);
        let poly = Polygon::new(ls, vec![]);
        // A point on the octagon
        let p = Point::new(5.0, 1.0);
        let dist = Euclidean.distance(&p, &poly);
        assert_relative_eq!(dist, 0.0);
    }
    #[test]
    // Point to Polygon, on boundary
    fn point_polygon_boundary_test2() {
        let exterior = LineString::from(vec![
            (0., 0.),
            (0., 0.0004),
            (0.0004, 0.0004),
            (0.0004, 0.),
            (0., 0.),
        ]);

        let poly = Polygon::new(exterior, vec![]);
        let bugged_point = Point::new(0.0001, 0.);
        assert_relative_eq!(Euclidean.distance(&poly, &bugged_point), 0.);
    }
    #[test]
    // Point to Polygon, empty Polygon
    fn point_polygon_empty_test() {
        // an empty Polygon
        let points = vec![];
        let ls = LineString::new(points);
        let poly = Polygon::new(ls, vec![]);
        // A point on the octagon
        let p = Point::new(2.5, 0.5);
        let dist = Euclidean.distance(&p, &poly);
        assert_relative_eq!(dist, 0.0);
    }
    #[test]
    // Point to Polygon with an interior ring
    fn point_polygon_interior_cutout_test() {
        // an octagon
        let ext_points = vec![
            (4., 1.),
            (5., 2.),
            (5., 3.),
            (4., 4.),
            (3., 4.),
            (2., 3.),
            (2., 2.),
            (3., 1.),
            (4., 1.),
        ];
        // cut out a triangle inside octagon
        let int_points = vec![(3.5, 3.5), (4.4, 1.5), (2.6, 1.5), (3.5, 3.5)];
        let ls_ext = LineString::from(ext_points);
        let ls_int = LineString::from(int_points);
        let poly = Polygon::new(ls_ext, vec![ls_int]);
        // A point inside the cutout triangle
        let p = Point::new(3.5, 2.5);
        let dist = Euclidean.distance(&p, &poly);

        // 0.41036467732879783 <-- Shapely
        assert_relative_eq!(dist, 0.41036467732879767);
    }

    #[test]
    fn line_distance_multipolygon_do_not_intersect_test() {
        // checks that the distance from the multipolygon
        // is equal to the distance from the closest polygon
        // taken in isolation, whatever that distance is
        let ls1 = LineString::from(vec![
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (5.0, 15.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ]);
        let ls2 = LineString::from(vec![
            (0.0, 30.0),
            (0.0, 25.0),
            (10.0, 25.0),
            (10.0, 30.0),
            (0.0, 30.0),
        ]);
        let ls3 = LineString::from(vec![
            (15.0, 30.0),
            (15.0, 25.0),
            (20.0, 25.0),
            (20.0, 30.0),
            (15.0, 30.0),
        ]);
        let pol1 = Polygon::new(ls1, vec![]);
        let pol2 = Polygon::new(ls2, vec![]);
        let pol3 = Polygon::new(ls3, vec![]);
        let mp = MultiPolygon::new(vec![pol1.clone(), pol2, pol3]);
        let pnt1 = Point::new(0.0, 15.0);
        let pnt2 = Point::new(10.0, 20.0);
        let ln = Line::new(pnt1.0, pnt2.0);
        let dist_mp_ln = Euclidean.distance(&ln, &mp);
        let dist_pol1_ln = Euclidean.distance(&ln, &pol1);
        assert_relative_eq!(dist_mp_ln, dist_pol1_ln);
    }

    #[test]
    fn point_distance_multipolygon_test() {
        let ls1 = LineString::from(vec![(0.0, 0.0), (1.0, 10.0), (2.0, 0.0), (0.0, 0.0)]);
        let ls2 = LineString::from(vec![(3.0, 0.0), (4.0, 10.0), (5.0, 0.0), (3.0, 0.0)]);
        let p1 = Polygon::new(ls1, vec![]);
        let p2 = Polygon::new(ls2, vec![]);
        let mp = MultiPolygon::new(vec![p1, p2]);
        let p = Point::new(50.0, 50.0);
        assert_relative_eq!(Euclidean.distance(&p, &mp), 60.959002616512684);
    }
    #[test]
    // Point to LineString
    fn point_linestring_distance_test() {
        // like an octagon, but missing the lowest horizontal segment
        let points = vec![
            (5., 1.),
            (4., 2.),
            (4., 3.),
            (5., 4.),
            (6., 4.),
            (7., 3.),
            (7., 2.),
            (6., 1.),
        ];
        let ls = LineString::from(points);
        // A Random point "inside" the LineString
        let p = Point::new(5.5, 2.1);
        let dist = Euclidean.distance(&p, &ls);
        assert_relative_eq!(dist, 1.1313708498984762);
    }
    #[test]
    // Point to LineString, point lies on the LineString
    fn point_linestring_contains_test() {
        // like an octagon, but missing the lowest horizontal segment
        let points = vec![
            (5., 1.),
            (4., 2.),
            (4., 3.),
            (5., 4.),
            (6., 4.),
            (7., 3.),
            (7., 2.),
            (6., 1.),
        ];
        let ls = LineString::from(points);
        // A point which lies on the LineString
        let p = Point::new(5.0, 4.0);
        let dist = Euclidean.distance(&p, &ls);
        assert_relative_eq!(dist, 0.0);
    }
    #[test]
    // Point to LineString, closed triangle
    fn point_linestring_triangle_test() {
        let points = vec![(3.5, 3.5), (4.4, 2.0), (2.6, 2.0), (3.5, 3.5)];
        let ls = LineString::from(points);
        let p = Point::new(3.5, 2.5);
        let dist = Euclidean.distance(&p, &ls);
        assert_relative_eq!(dist, 0.5);
    }
    #[test]
    // Point to LineString, empty LineString
    fn point_linestring_empty_test() {
        let points = vec![];
        let ls = LineString::new(points);
        let p = Point::new(5.0, 4.0);
        let dist = Euclidean.distance(&p, &ls);
        assert_relative_eq!(dist, 0.0);
    }
    #[test]
    fn distance_multilinestring_test() {
        let v1 = LineString::from(vec![(0.0, 0.0), (1.0, 10.0)]);
        let v2 = LineString::from(vec![(1.0, 10.0), (2.0, 0.0), (3.0, 1.0)]);
        let mls = MultiLineString::new(vec![v1, v2]);
        let p = Point::new(50.0, 50.0);
        assert_relative_eq!(Euclidean.distance(&p, &mls), 63.25345840347388);
    }
    #[test]
    fn distance1_test() {
        assert_relative_eq!(
            Euclidean.distance(&Point::new(0., 0.), &Point::new(1., 0.)),
            1.
        );
    }
    #[test]
    fn distance2_test() {
        let dist = Euclidean.distance(&Point::new(-72.1235, 42.3521), &Point::new(72.1260, 70.612));
        assert_relative_eq!(dist, 146.99163308930207);
    }
    #[test]
    fn distance_multipoint_test() {
        let v = vec![
            Point::new(0.0, 10.0),
            Point::new(1.0, 1.0),
            Point::new(10.0, 0.0),
            Point::new(1.0, -1.0),
            Point::new(0.0, -10.0),
            Point::new(-1.0, -1.0),
            Point::new(-10.0, 0.0),
            Point::new(-1.0, 1.0),
            Point::new(0.0, 10.0),
        ];
        let mp = MultiPoint::new(v);
        let p = Point::new(50.0, 50.0);
        assert_relative_eq!(Euclidean.distance(&p, &mp), 64.03124237432849)
    }
    #[test]
    fn distance_line_test() {
        let line0 = Line::from([(0., 0.), (5., 0.)]);
        let p0 = Point::new(2., 3.);
        let p1 = Point::new(3., 0.);
        let p2 = Point::new(6., 0.);
        assert_relative_eq!(Euclidean.distance(&line0, &p0), 3.);
        assert_relative_eq!(Euclidean.distance(&p0, &line0), 3.);

        assert_relative_eq!(Euclidean.distance(&line0, &p1), 0.);
        assert_relative_eq!(Euclidean.distance(&p1, &line0), 0.);

        assert_relative_eq!(Euclidean.distance(&line0, &p2), 1.);
        assert_relative_eq!(Euclidean.distance(&p2, &line0), 1.);
    }
    #[test]
    fn distance_line_line_test() {
        let line0 = Line::from([(0., 0.), (5., 0.)]);
        let line1 = Line::from([(2., 1.), (7., 2.)]);
        assert_relative_eq!(Euclidean.distance(&line0, &line1), 1.);
        assert_relative_eq!(Euclidean.distance(&line1, &line0), 1.);
    }
    #[test]
    // See https://github.com/georust/geo/issues/476
    fn distance_line_polygon_test() {
        let line = Line::new(
            coord! {
                x: -0.17084137691985102,
                y: 0.8748085493016657,
            },
            coord! {
                x: -0.17084137691985102,
                y: 0.09858870312437906,
            },
        );
        let poly: Polygon<f64> = polygon![
            coord! {
                x: -0.10781391405721802,
                y: -0.15433610862574643,
            },
            coord! {
                x: -0.7855276236615211,
                y: 0.23694208404779793,
            },
            coord! {
                x: -0.7855276236615214,
                y: -0.5456143012992907,
            },
            coord! {
                x: -0.10781391405721802,
                y: -0.15433610862574643,
            },
        ];
        assert_eq!(Euclidean.distance(&line, &poly), 0.18752558079168907);
    }
    #[test]
    // test edge-vertex minimum distance
    fn test_minimum_polygon_distance() {
        let points_raw = [
            (126., 232.),
            (126., 212.),
            (112., 202.),
            (97., 204.),
            (87., 215.),
            (87., 232.),
            (100., 246.),
            (118., 247.),
        ];
        let points = points_raw
            .iter()
            .map(|e| Point::new(e.0, e.1))
            .collect::<Vec<_>>();
        let poly1 = Polygon::new(LineString::from(points), vec![]);

        let points_raw_2 = [
            (188., 231.),
            (189., 207.),
            (174., 196.),
            (164., 196.),
            (147., 220.),
            (158., 242.),
            (177., 242.),
        ];
        let points2 = points_raw_2
            .iter()
            .map(|e| Point::new(e.0, e.1))
            .collect::<Vec<_>>();
        let poly2 = Polygon::new(LineString::from(points2), vec![]);
        let dist = nearest_neighbour_distance(poly1.exterior(), poly2.exterior());
        assert_relative_eq!(dist, 21.0);
    }
    #[test]
    // test vertex-vertex minimum distance
    fn test_minimum_polygon_distance_2() {
        let points_raw = [
            (118., 200.),
            (153., 179.),
            (106., 155.),
            (88., 190.),
            (118., 200.),
        ];
        let points = points_raw
            .iter()
            .map(|e| Point::new(e.0, e.1))
            .collect::<Vec<_>>();
        let poly1 = Polygon::new(LineString::from(points), vec![]);

        let points_raw_2 = [
            (242., 186.),
            (260., 146.),
            (182., 175.),
            (216., 193.),
            (242., 186.),
        ];
        let points2 = points_raw_2
            .iter()
            .map(|e| Point::new(e.0, e.1))
            .collect::<Vec<_>>();
        let poly2 = Polygon::new(LineString::from(points2), vec![]);
        let dist = nearest_neighbour_distance(poly1.exterior(), poly2.exterior());
        assert_relative_eq!(dist, 29.274562336608895);
    }
    #[test]
    // test edge-edge minimum distance
    fn test_minimum_polygon_distance_3() {
        let points_raw = [
            (182., 182.),
            (182., 168.),
            (138., 160.),
            (136., 193.),
            (182., 182.),
        ];
        let points = points_raw
            .iter()
            .map(|e| Point::new(e.0, e.1))
            .collect::<Vec<_>>();
        let poly1 = Polygon::new(LineString::from(points), vec![]);

        let points_raw_2 = [
            (232., 196.),
            (234., 150.),
            (194., 165.),
            (194., 191.),
            (232., 196.),
        ];
        let points2 = points_raw_2
            .iter()
            .map(|e| Point::new(e.0, e.1))
            .collect::<Vec<_>>();
        let poly2 = Polygon::new(LineString::from(points2), vec![]);
        let dist = nearest_neighbour_distance(poly1.exterior(), poly2.exterior());
        assert_relative_eq!(dist, 12.0);
    }
    #[test]
    fn test_large_polygon_distance() {
        let ls = geo_test_fixtures::norway_main::<f64>();
        let poly1 = Polygon::new(ls, vec![]);
        let vec2 = vec![
            (4.921875, 66.33750501996518),
            (3.69140625, 65.21989393613207),
            (6.15234375, 65.07213008560697),
            (4.921875, 66.33750501996518),
        ];
        let poly2 = Polygon::new(vec2.into(), vec![]);
        let distance = Euclidean.distance(&poly1, &poly2);
        // GEOS says 2.2864896295566055
        assert_relative_eq!(distance, 2.2864896295566055);
    }
    #[test]
    // A polygon inside another polygon's ring; they're disjoint in the DE-9IM sense:
    // FF2FF1212
    fn test_poly_in_ring() {
        let shell = geo_test_fixtures::shell::<f64>();
        let ring = geo_test_fixtures::ring::<f64>();
        let poly_in_ring = geo_test_fixtures::poly_in_ring::<f64>();
        // inside is "inside" outside's ring, but they are disjoint
        let outside = Polygon::new(shell, vec![ring]);
        let inside = Polygon::new(poly_in_ring, vec![]);
        assert_relative_eq!(Euclidean.distance(&outside, &inside), 5.992772737231033);
    }
    #[test]
    // two ring LineStrings; one encloses the other but they neither touch nor intersect
    fn test_linestring_distance() {
        let ring = geo_test_fixtures::ring::<f64>();
        let poly_in_ring = geo_test_fixtures::poly_in_ring::<f64>();
        assert_relative_eq!(Euclidean.distance(&ring, &poly_in_ring), 5.992772737231033);
    }
    #[test]
    // Line-Polygon test: closest point on Polygon is NOT nearest to a Line end-point
    fn test_line_polygon_simple() {
        let line = Line::from([(0.0, 0.0), (0.0, 3.0)]);
        let v = vec![(5.0, 1.0), (5.0, 2.0), (0.25, 1.5), (5.0, 1.0)];
        let poly = Polygon::new(v.into(), vec![]);
        assert_relative_eq!(Euclidean.distance(&line, &poly), 0.25);
    }
    #[test]
    // Line-Polygon test: Line intersects Polygon
    fn test_line_polygon_intersects() {
        let line = Line::from([(0.5, 0.0), (0.0, 3.0)]);
        let v = vec![(5.0, 1.0), (5.0, 2.0), (0.25, 1.5), (5.0, 1.0)];
        let poly = Polygon::new(v.into(), vec![]);
        assert_relative_eq!(Euclidean.distance(&line, &poly), 0.0);
    }
    #[test]
    // Line-Polygon test: Line contained by interior ring
    fn test_line_polygon_inside_ring() {
        let line = Line::from([(4.4, 1.5), (4.45, 1.5)]);
        let v = vec![(5.0, 1.0), (5.0, 2.0), (0.25, 1.0), (5.0, 1.0)];
        let v2 = vec![(4.5, 1.2), (4.5, 1.8), (3.5, 1.2), (4.5, 1.2)];
        let poly = Polygon::new(v.into(), vec![v2.into()]);
        assert_relative_eq!(Euclidean.distance(&line, &poly), 0.04999999999999982);
    }
    #[test]
    // LineString-Line test
    fn test_linestring_line_distance() {
        let line = Line::from([(0.0, 0.0), (0.0, 2.0)]);
        let ls: LineString<_> = vec![(3.0, 0.0), (1.0, 1.0), (3.0, 2.0)].into();
        assert_relative_eq!(Euclidean.distance(&ls, &line), 1.0);
    }

    #[test]
    // Triangle-Point test: point on vertex
    fn test_triangle_point_on_vertex_distance() {
        let triangle = Triangle::from([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)]);
        let point = Point::new(0.0, 0.0);
        assert_relative_eq!(Euclidean.distance(&triangle, &point), 0.0);
    }

    #[test]
    // Triangle-Point test: point on edge
    fn test_triangle_point_on_edge_distance() {
        let triangle = Triangle::from([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)]);
        let point = Point::new(1.5, 0.0);
        assert_relative_eq!(Euclidean.distance(&triangle, &point), 0.0);
    }

    #[test]
    // Triangle-Point test
    fn test_triangle_point_distance() {
        let triangle = Triangle::from([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)]);
        let point = Point::new(2.0, 3.0);
        assert_relative_eq!(Euclidean.distance(&triangle, &point), 1.0);
    }

    #[test]
    // Triangle-Point test: point within triangle
    fn test_triangle_point_inside_distance() {
        let triangle = Triangle::from([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)]);
        let point = Point::new(1.0, 0.5);
        assert_relative_eq!(Euclidean.distance(&triangle, &point), 0.0);
    }

    #[test]
    fn convex_and_nearest_neighbour_comparison() {
        let ls1: LineString<f64> = vec![
            Coord::from((57.39453770777941, 307.60533608924663)),
            Coord::from((67.1800355576469, 309.6654408997451)),
            Coord::from((84.89693692793338, 225.5101593908847)),
            Coord::from((75.1114390780659, 223.45005458038628)),
            Coord::from((57.39453770777941, 307.60533608924663)),
        ]
        .into();
        let first_polygon: Polygon<f64> = Polygon::new(ls1, vec![]);
        let ls2: LineString<f64> = vec![
            Coord::from((138.11769866645008, -45.75134112915392)),
            Coord::from((130.50230476949187, -39.270154833870336)),
            Coord::from((184.94426964987397, 24.699153900578573)),
            Coord::from((192.55966354683218, 18.217967605294987)),
            Coord::from((138.11769866645008, -45.75134112915392)),
        ]
        .into();
        let second_polygon = Polygon::new(ls2, vec![]);

        assert_relative_eq!(
            Euclidean.distance(&first_polygon, &second_polygon),
            224.35357967013238
        );
    }
    #[test]
    fn fast_path_regression() {
        // this test will fail if the fast path algorithm is reintroduced without being fixed
        let p1 = polygon!(
            (x: 0_f64, y: 0_f64),
            (x: 300_f64, y: 0_f64),
            (x: 300_f64, y: 100_f64),
            (x: 0_f64, y: 100_f64),
        )
        .orient(Direction::Default);
        let p2 = polygon!(
            (x: 100_f64, y: 150_f64),
            (x: 150_f64, y: 200_f64),
            (x: 50_f64, y: 200_f64),
        )
        .orient(Direction::Default);
        let p3 = polygon!(
            (x: 0_f64, y: 0_f64),
            (x: 300_f64, y: 0_f64),
            (x: 300_f64, y: 100_f64),
            (x: 0_f64, y: 100_f64),
        )
        .orient(Direction::Reversed);
        let p4 = polygon!(
            (x: 100_f64, y: 150_f64),
            (x: 150_f64, y: 200_f64),
            (x: 50_f64, y: 200_f64),
        )
        .orient(Direction::Reversed);
        assert_eq!(Euclidean.distance(&p1, &p2), 50.0f64);
        assert_eq!(Euclidean.distance(&p3, &p4), 50.0f64);
        assert_eq!(Euclidean.distance(&p1, &p4), 50.0f64);
        assert_eq!(Euclidean.distance(&p2, &p3), 50.0f64);
    }
    #[test]
    fn all_types_geometry_collection_test() {
        let p = Point::new(0.0, 0.0);
        let line = Line::from([(-1.0, -1.0), (-2.0, -2.0)]);
        let ls = LineString::from(vec![(0.0, 0.0), (1.0, 10.0), (2.0, 0.0)]);
        let poly = Polygon::new(
            LineString::from(vec![(0.0, 0.0), (1.0, 10.0), (2.0, 0.0), (0.0, 0.0)]),
            vec![],
        );
        let tri = Triangle::from([(0.0, 0.0), (1.0, 10.0), (2.0, 0.0)]);
        let rect = Rect::new((0.0, 0.0), (-1.0, -1.0));

        let ls1 = LineString::from(vec![(0.0, 0.0), (1.0, 10.0), (2.0, 0.0), (0.0, 0.0)]);
        let ls2 = LineString::from(vec![(3.0, 0.0), (4.0, 10.0), (5.0, 0.0), (3.0, 0.0)]);
        let p1 = Polygon::new(ls1, vec![]);
        let p2 = Polygon::new(ls2, vec![]);
        let mpoly = MultiPolygon::new(vec![p1, p2]);

        let v = vec![
            Point::new(0.0, 10.0),
            Point::new(1.0, 1.0),
            Point::new(10.0, 0.0),
            Point::new(1.0, -1.0),
            Point::new(0.0, -10.0),
            Point::new(-1.0, -1.0),
            Point::new(-10.0, 0.0),
            Point::new(-1.0, 1.0),
            Point::new(0.0, 10.0),
        ];
        let mpoint = MultiPoint::new(v);

        let v1 = LineString::from(vec![(0.0, 0.0), (1.0, 10.0)]);
        let v2 = LineString::from(vec![(1.0, 10.0), (2.0, 0.0), (3.0, 1.0)]);
        let mls = MultiLineString::new(vec![v1, v2]);

        let gc = GeometryCollection(vec![
            Geometry::Point(p),
            Geometry::Line(line),
            Geometry::LineString(ls),
            Geometry::Polygon(poly),
            Geometry::MultiPoint(mpoint),
            Geometry::MultiLineString(mls),
            Geometry::MultiPolygon(mpoly),
            Geometry::Triangle(tri),
            Geometry::Rect(rect),
        ]);

        let test_p = Point::new(50., 50.);
        assert_relative_eq!(Euclidean.distance(&test_p, &gc), 60.959002616512684);

        let test_multipoint = MultiPoint::new(vec![test_p]);
        assert_relative_eq!(
            Euclidean.distance(&test_multipoint, &gc),
            60.959002616512684
        );

        let test_line = Line::from([(50., 50.), (60., 60.)]);
        assert_relative_eq!(Euclidean.distance(&test_line, &gc), 60.959002616512684);

        let test_ls = LineString::from(vec![(50., 50.), (60., 60.), (70., 70.)]);
        assert_relative_eq!(Euclidean.distance(&test_ls, &gc), 60.959002616512684);

        let test_mls = MultiLineString::new(vec![test_ls]);
        assert_relative_eq!(Euclidean.distance(&test_mls, &gc), 60.959002616512684);

        let test_poly = Polygon::new(
            LineString::from(vec![
                (50., 50.),
                (60., 50.),
                (60., 60.),
                (55., 55.),
                (50., 50.),
            ]),
            vec![],
        );
        assert_relative_eq!(Euclidean.distance(&test_poly, &gc), 60.959002616512684);

        let test_multipoly = MultiPolygon::new(vec![test_poly]);
        assert_relative_eq!(Euclidean.distance(&test_multipoly, &gc), 60.959002616512684);

        let test_tri = Triangle::from([(50., 50.), (60., 50.), (55., 55.)]);
        assert_relative_eq!(Euclidean.distance(&test_tri, &gc), 60.959002616512684);

        let test_rect = Rect::new(coord! { x: 50., y: 50. }, coord! { x: 60., y: 60. });
        assert_relative_eq!(Euclidean.distance(&test_rect, &gc), 60.959002616512684);

        let test_gc = GeometryCollection(vec![Geometry::Rect(test_rect)]);
        assert_relative_eq!(Euclidean.distance(&test_gc, &gc), 60.959002616512684);
    }
}

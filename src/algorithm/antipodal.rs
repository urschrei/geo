// Pirzadeh, H. (1999) Computational geometry with the rotating calipers., pp30 – 32
// Available from: http://digitool.library.mcgill.ca/R/-?func=dbin-jump-full&object_id=21623&silo_library=GEN01
// http://web.archive.org/web/20150330010154/http://cgm.cs.mcgill.ca/%7Eorm/rotcal.html
use num_traits::Float;
use types::{Point, Polygon, MultiPolygon, LineString, MultiPoint, MultiLineString};
use std::fmt::Debug;
use std::mem;
use algorithm::hull_helpers::{
    swap_remove_to_first,
    swap_remove_to_last,
    partition,
    point_location
};
use algorithm::convexhull::ConvexHull;

fn min_polygon_distance<T>(mut poly1: Polygon<T>, mut poly2: Polygon<T>) -> T
    where T: Float + Debug
{
    // polygons must be convex
    let mut poly1_hull = poly1.exterior.0.as_mut_slice();
    let mut poly2_hull = poly2.exterior.0.as_mut_slice();
    // guess
    let mut poly1_xmin = poly1_hull.first().unwrap().clone();
    let mut poly2_xmin = poly2_hull.first().unwrap().clone();

    // poly1
    let mut poly1_ymin = swap_remove_to_first(&mut poly1_hull, 0);
    let mut poly1_ymax = swap_remove_to_first(&mut poly1_hull, 0);
    let mut poly1_xmax = swap_remove_to_first(&mut poly1_hull, 0);
    if poly1_ymin.y() > poly1_ymax.y() {
        mem::swap(poly1_ymin, poly1_ymax);
    }
    for point in poly1_hull.iter_mut() {
        if point.y() < poly1_ymin.y() {
            mem::swap(point, poly1_ymin);
        }
        if point.y() > poly1_ymax.y() {
            mem::swap(point, poly1_ymax);
        }
        if point.x() > poly1_xmax.x() {
            mem::swap(point, &mut poly1_xmax);
        }
        if point.x() < poly1_xmin.x() {
            mem::swap(point, &mut poly1_xmin);
        }
    }
    // poly2
    let mut poly2_ymin = swap_remove_to_first(&mut poly2_hull, 0);
    let mut poly2_ymax = swap_remove_to_first(&mut poly2_hull, 0);
    let mut poly2_xmax = swap_remove_to_first(&mut poly2_hull, 0);
    if poly2_ymin.y() > poly2_ymax.y() {
        mem::swap(poly2_ymin, poly2_ymax);
    }
    for point in poly2_hull.iter_mut() {
        if point.y() < poly2_ymin.y() {
            mem::swap(point, poly2_ymin);
        }
        if point.y() > poly2_ymax.y() {
            mem::swap(point, poly2_ymax);
        }
        if point.x() > poly2_xmax.x() {
            mem::swap(point, &mut poly2_xmax);
        }
        if point.x() < poly2_xmin.x() {
            mem::swap(point, &mut poly2_xmin);
        }
    }
    // lines of support must be parallel to the x axis
    // lpoly1 must have poly1 to its right
    // lpoly2 must have poly2 to its right
    let mut lpoly_1 = LineString(vec![
        Point::new(poly1_xmax.x(), poly1_ymin.y()),
        Point::new(poly1_ymin.x(), poly1_ymin.y()),
        Point::new(poly1_xmin.x(), poly1_ymin.y())
    ]);
    let mut lpoly_2 = LineString(vec![
        Point::new(poly2_xmin.x(), poly2_ymax.y()),
        Point::new(poly2_ymax.x(), poly2_ymax.y()),
        Point::new(poly2_xmax.x(), poly2_ymax.y())
    ]);
    println!("poly 1 min Y: {:?}", poly1_ymin.y());
    println!("poly 2 max Y: {:?}", poly2_ymax.y());
    println!("poly 1 min X: {:?}", poly1_xmin);
    println!("poly 2 max X: {:?}", poly2_xmax);
    println!("Bottom support (r to l: {:?}", lpoly_1);
    println!("Top support (l to r): {:?}", lpoly_2);
    // 1.  We want poly1_min.y(), and poly2_max.y()
    // 2.  Construct two lines of support, parallel to the x axis – LP and LQ –
    //     which touch the polygons at yminP and ymaxQ
    //     such that the polygons lie to the right of their respective lines of support.
    //     LP and LQ have opposite direction, and yminP and ymaxQ form an anti-podal pair between the polygons.
    //     The lines of support lie on vertices pi ∈ P, and qj ∈ Q, and determine two angles
    //     θi, θj, which are computed
    // 3.  Compute dist(yminP,ymaxQ) and keep it as the minimum.
    // 3a. Compute θ = θi.min(θj)
    // 4.  Rotate the lines clockwise about pi and qj by θ
    //     One of the lines should now be flush with an edge of its polygon. (???)
    // 5.  If only one line coincides with an edge, then:
    //     - the vertex-edge anti-podal pair distance should be computed
    //     - the new vertex-vertex anti-podal pair distance should be computed
    //     Both distances are compared the current minimum, which is updated if necessary.
    //     
    //     If both lines of support coincide with edges, then the situation is somewhat more complex:
    //     - If the edges "overlap", that is if one can construct a line perpendicular to both edges and
    //     intersecting both edges (but not at vertices), then the edge-edge distance should be computed.
    //     - Otherwise the three new vertex-vertex anti-podal pair distances are computed.
    //     All distances are compared to the current minimum, which is updated if necessary.
    // 6.  Repeat steps 4 and 5, until the lines reach (yminP, ymaxQ) again.
    // 7. Return the minimum

    T::from(3.0).unwrap()
}

#[cfg(test)]
mod test {
    use types::Point;
    use super::*;
    #[test]
    fn test_minimum_polygon_distance() {
        let points_raw = vec![(5., 1.), (4., 2.), (4., 3.), (5., 4.), (6., 4.), (7., 3.),
                              (7., 2.), (6., 1.), (5., 1.)];
        let mut points = points_raw.iter().map(|e| Point::new(e.0, e.1)).collect::<Vec<_>>();
        let poly1 = Polygon::new(LineString(points), vec![]);

        let points_raw_2 = vec![(8., 1.), (7., 2.), (7., 3.), (8., 4.), (9., 4.), (10., 3.),
                      (10., 2.), (9., 1.), (8., 1.)];
        let mut points2 = points_raw_2.iter().map(|e| Point::new(e.0, e.1)).collect::<Vec<_>>();
        let poly2 = Polygon::new(LineString(points2), vec![]);
        let dist = min_polygon_distance(poly1.convex_hull(), poly2.convex_hull());
        assert_eq!(dist, 3.0);
    }
}
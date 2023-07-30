use crate::{
    relate::geomgraph::Quadrant, Coord, CoordNum, CoordsIter, GeoFloat, GeoNum, HasKernel, Point,
};

// Ported from JTS on 30/07/2023
// path: jts/modules/core/src/main/java/org/locationtech/jts/geomgraph/index/MonotoneChainIndexer.java
// latest JTS repo file commit at time of port: 5d22d87

/// Monotone chains are a way of partitioning the segments of an edge to
/// allow for fast searching of intersections.
/// Specifically, a sequence of contiguous line segments
/// is a monotone chain if all the vectors defined by the oriented segments
/// lies in the same quadrant.
///
/// Monotone Chains have the following useful properties:
///
/// - the segments within a monotone chain will never intersect each other
///  - the envelope of any contiguous subset of the segments in a monotone chain is simply the envelope of the endpoints of the subset.
///
/// Property 1 means that there is no need to test pairs of segments from within
/// the same monotone chain for intersection.
///
/// Property 2 allows
/// binary search to be used to find the intersection points of two monotone chains.
/// For many types of real-world data, these properties eliminate a large number of
/// segment comparisons, producing substantial speed gains.
///
/// Note that due to the efficient intersection test, there is no need to limit the size
/// of chains to obtain fast performance.
pub(crate) struct MonotoneChainIndexer {}

impl MonotoneChainIndexer {
    /// This is simple in Rust, but is included here for completeness
    pub(crate) fn to_int_array(it: impl Iterator) -> Vec<usize> {
        it.enumerate().map(|(idx, _)| idx).collect()
    }

    /// find the startpoint (and endpoints) of all monotone chains in this edge
    pub(crate) fn chain_start_indices<T>(self, container: &[Coord<T>]) -> Vec<usize>
    where
        T: GeoFloat,
    {
        let start = 0;
        let mut start_index_vec = Vec::with_capacity(container.len() / 2);
        // use heuristic to size initial array
        start_index_vec.push(start);
        while start < container.len() {
            let last = Self::find_chain_end(container, start);
            start_index_vec.push(last);
        }
        start_index_vec
    }

    fn find_chain_end<T>(container: &[Coord<T>], start: usize) -> usize
    where
        T: GeoFloat,
    {
        let first = container[start];
        let mut last = start + 1;
        let end = container[last];

        // TODO: handle None
        let chain_quad = Quadrant::new_from_coords(first, end).unwrap();
        while last < container.len() {
            // TODO: handle None
            let quad = Quadrant::new_from_coords(container[last - 1], container[last]).unwrap();
            if quad != chain_quad {
                break;
            }
            last -= 1;
        }
        last - 1
    }
}

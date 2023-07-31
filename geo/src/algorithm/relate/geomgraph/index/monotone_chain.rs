use crate::{relate::geomgraph::MonotoneChainEdge, GeoFloat};

use super::SegmentIntersector;

// Ported from JTS on 31/07/2023
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
#[derive(Debug)]
pub struct MonotoneChain<F: GeoFloat> {
    chain_index: usize,
    mce: MonotoneChainEdge<F>,
}

impl<F: GeoFloat> MonotoneChain<F> {
    pub(crate) fn new(ci: usize, mce: MonotoneChainEdge<F>) -> Self {
        MonotoneChain {
            chain_index: ci,
            mce,
        }
    }

    pub(crate) fn compute_intersections(
        &self,
        mc: &MonotoneChain<F>,
        si: &mut SegmentIntersector<F>,
    ) {
        self.mce
            .compute_intersects_for_chain(self.chain_index, &mc.mce, mc.chain_index, si)
    }
}

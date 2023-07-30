use std::cell::RefCell;

use rstar::Envelope;

use super::{index::MonotoneChainIndexer, index::SegmentIntersector, Edge};
use crate::{Coord, GeoFloat};

// Ported from JTS on 30/07/2023
// path: jts/modules/core/src/main/java/org/locationtech/jts/geomgraph/index/MonotoneChainIndexer.java
// latest JTS repo file commit at time of port: 5d22d87

/// A one-dimensional line in a Monotone Chain
#[derive(Debug)]
pub(crate) struct MonotoneChainEdge<F: GeoFloat> {
    /// JTS calls this `pts` and says it's a "cached reference, for efficiency".
    ///
    /// I have no idea what that means. Also consider that `Edge` only has 2 Coords anyway and should
    /// perhaps be a Copy type ?
    coords: Vec<Coord<F>>,
    edge: RefCell<Edge<F>>,
    /// the lists of start/end indexes of the monotone chains.
    /// Includes the end point of the edge as a sentinel
    startindex: Vec<usize>,
}

impl<F: GeoFloat> MonotoneChainEdge<F> {
    pub(crate) fn new(e: Edge<F>) -> Self {
        let mcb = MonotoneChainIndexer {};
        let coords = e.coords().to_vec();
        let si = mcb.chain_start_indices(&coords);
        MonotoneChainEdge {
            coords,
            edge: RefCell::new(e),
            startindex: si,
        }
    }

    pub(crate) fn coords(&self) -> &[Coord<F>] {
        &self.coords
    }

    pub(crate) fn start_indexes(&self) -> &[usize] {
        &self.startindex
    }

    pub(crate) fn min_x(&self, chain_index: usize) -> F {
        let x1 = self.coords[self.startindex[chain_index]].x;
        let x2 = self.coords[self.startindex[chain_index + 1]].x;
        x1.min(x2)
    }

    pub(crate) fn max_x(&self, chain_index: usize) -> F {
        let x1 = self.coords[self.startindex[chain_index]].x;
        let x2 = self.coords[self.startindex[chain_index + 1]].x;
        x1.max(x2)
    }

    /// Compute all intersections for a given Monotone Chain
    // It's not clear to me why JTS has both a public and private function with this name,
    // with the private function accepting two additional arguments (start0 and end0)
    // I've renamed it compute_actual here
    pub(crate) fn compute_intersects_for_chain(
        &self,
        chain_index0: usize,
        mce: &MonotoneChainEdge<F>,
        chain_index1: usize,
        si: &mut SegmentIntersector<F>,
    ) {
        let s1 = mce.startindex[chain_index1];
        let e1 = mce.startindex[chain_index1 + 1];
        self.compute_actual(
            self.startindex[chain_index0],
            self.startindex[chain_index0 + 1],
            mce,
            s1,
            e1,
            si,
        );
    }

    fn compute_actual(
        &self,
        start0: usize,
        end0: usize,
        mce: &MonotoneChainEdge<F>,
        start1: usize,
        end1: usize,
        ei: &mut SegmentIntersector<F>,
    ) {
        // terminating condition for the recursion
        if end0 - start0 == 1 && end1 - start1 == 1 {
            // ei.add_intersections(edge0, segment_index_0, edge1, segment_index_1)
            ei.add_intersections(&self.edge, start0, &mce.edge, start1);
            return;
        }
        // nothing to do if the envelopes of these chains don't overlap
        if !self.overlaps(start0, end0, mce, start1, end1) {
            return;
        }
        // the chains overlap, so split each in half and iterate  (binary search)
        let mid0 = (start0 + end0) / 2;
        let mid1 = (start1 + end1) / 2;
        if start0 < mid0 {
            if start1 < mid1 {
                self.compute_actual(start0, mid0, mce, start1, mid1, ei);
            }
            if mid1 < end1 {
                self.compute_actual(start0, mid0, mce, mid1, end1, ei);
            }
        }
    }
    /// Tests whether the envelopes of two chain sections overlap (intersect).
    fn overlaps(
        &self,
        start0: usize,
        end0: usize,
        mce: &MonotoneChainEdge<F>,
        start1: usize,
        end1: usize,
    ) -> bool {
        let env_a = rstar::AABB::from_corners(self.coords[start0], self.coords[end0]);
        let env_b = rstar::AABB::from_corners(mce.coords[start1], mce.coords[end1]);
        env_a.intersects(&env_b)
    }
}

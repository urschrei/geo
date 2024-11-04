mod i_overlay_integration;
#[cfg(test)]
mod tests;

pub use i_overlay_integration::BoolOpsNum;

use crate::geometry::{LineString, MultiLineString, MultiPolygon, Polygon};
use rstar::{primitives::CachedEnvelope, ParentNode, RTree, RTreeNode, RTreeObject};

/// Boolean Operations on geometry.
///
/// Boolean operations are set operations on geometries considered as a subset
/// of the 2-D plane. The operations supported are: intersection, union, xor or
/// symmetric difference, and set-difference on pairs of 2-D geometries and
/// clipping a 1-D geometry with self.
///
/// These operations are implemented on [`Polygon`] and the [`MultiPolygon`]
/// geometries.
///
/// # Validity
///
/// Note that the operations are strictly well-defined only on *valid*
/// geometries. However, the implementation generally works well as long as the
/// interiors of polygons are contained within their corresponding exteriors.
///
/// Degenerate 2-d geoms with 0 area are handled, and ignored by the algorithm.
/// In particular, taking `union` with an empty geom should remove degeneracies
/// and fix invalid polygons as long the interior-exterior requirement above is
/// satisfied.
///
/// # Performance
///
/// For union operations on a collection of overlapping and / or adjacent [`Polygon`]s
/// (e.g. contained in a `Vec` or a [`MultiPolygon`]), using [`UnaryUnion`] will
/// yield far better performance.
pub trait BooleanOps {
    type Scalar: BoolOpsNum;

    /// The exterior and interior rings of the geometry.
    ///
    /// It doesn't particularly matter which order they are in, as the topology algorithm counts crossings
    /// to determine the interior and exterior of the polygon.
    ///
    /// It is required that the rings are from valid geometries, that the rings not overlap.
    /// In the case of a MultiPolygon, this requires that none of its polygon's interiors may overlap.
    fn rings(&self) -> impl Iterator<Item = &LineString<Self::Scalar>>;

    fn boolean_op(
        &self,
        other: &impl BooleanOps<Scalar = Self::Scalar>,
        op: OpType,
    ) -> MultiPolygon<Self::Scalar> {
        use i_overlay::core::fill_rule::FillRule;
        use i_overlay::core::overlay::ShapeType;
        use i_overlay_integration::{convert, BoolOpsOverlay, BoolOpsOverlayGraph};
        let mut overlay = <Self::Scalar as BoolOpsNum>::OverlayType::new();

        for ring in self.rings() {
            overlay.add_path(convert::ring_to_shape_path(ring), ShapeType::Subject);
        }
        for ring in other.rings() {
            overlay.add_path(convert::ring_to_shape_path(ring), ShapeType::Clip);
        }

        let graph = overlay.into_graph(FillRule::EvenOdd);
        let shapes = graph.extract_shapes(op.into());

        convert::multi_polygon_from_shapes(shapes)
    }

    fn intersection(
        &self,
        other: &impl BooleanOps<Scalar = Self::Scalar>,
    ) -> MultiPolygon<Self::Scalar> {
        self.boolean_op(other, OpType::Intersection)
    }
    fn union(&self, other: &impl BooleanOps<Scalar = Self::Scalar>) -> MultiPolygon<Self::Scalar> {
        self.boolean_op(other, OpType::Union)
    }
    fn xor(&self, other: &impl BooleanOps<Scalar = Self::Scalar>) -> MultiPolygon<Self::Scalar> {
        self.boolean_op(other, OpType::Xor)
    }
    fn difference(
        &self,
        other: &impl BooleanOps<Scalar = Self::Scalar>,
    ) -> MultiPolygon<Self::Scalar> {
        self.boolean_op(other, OpType::Difference)
    }

    /// Clip a 1-D geometry with self.
    ///
    /// Returns the portion of `ls` that lies within `self` (known as the set-theoeretic
    /// intersection) if `invert` is false, and the difference (`ls - self`) otherwise.
    fn clip(
        &self,
        multi_line_string: &MultiLineString<Self::Scalar>,
        invert: bool,
    ) -> MultiLineString<Self::Scalar> {
        use i_overlay::core::fill_rule::FillRule;
        use i_overlay::string::clip::ClipRule;
        use i_overlay_integration::{convert, BoolOpsStringGraph, BoolOpsStringOverlay};

        let mut overlay = <Self::Scalar as BoolOpsNum>::StringOverlayType::new();

        for ring in self.rings() {
            overlay.add_shape_path(convert::ring_to_shape_path(ring));
        }
        for line_string in multi_line_string {
            for line in line_string.lines() {
                let line = [
                    Self::Scalar::to_bops_coord(line.start),
                    Self::Scalar::to_bops_coord(line.end),
                ];
                overlay.add_string_line(line)
            }
        }

        let graph = overlay.into_graph(FillRule::EvenOdd);
        let paths = graph.clip_string_lines(ClipRule {
            invert,
            boundary_included: true,
        });
        convert::multi_line_string_from_paths(paths)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum OpType {
    Intersection,
    Union,
    Difference,
    Xor,
}

/// Efficient [BooleanOps::union] of adjacent / overlapping geometries
///
/// For geometries with a high degree of overlap or adjacency
/// (for instance, merging a large contiguous area made up of many adjacent polygons)
/// this method will be orders of magnitude faster than a manual iteration and union approach.
pub trait UnaryUnion {
    type Scalar: BoolOpsNum;

    /// Construct a tree of all the input geometries and progressively union them from the "bottom up"
    ///
    /// This is considerably more efficient than using e.g. `fold()` over an iterator of Polygons.
    /// # Examples
    ///
    /// ```
    /// use geo::{BooleanOps, UnaryUnion};
    /// use geo::{MultiPolygon, polygon};
    /// let poly1 = polygon![
    ///     (x: 0.0, y: 0.0),
    ///     (x: 4.0, y: 0.0),
    ///     (x: 4.0, y: 4.0),
    ///     (x: 0.0, y: 4.0),
    ///     (x: 0.0, y: 0.0),
    /// ];
    /// let poly2 = polygon![
    ///     (x: 4.0, y: 0.0),
    ///     (x: 8.0, y: 0.0),
    ///     (x: 8.0, y: 4.0),
    ///     (x: 4.0, y: 4.0),
    ///     (x: 4.0, y: 0.0),
    /// ];
    /// let merged = &poly1.union(&poly2);
    /// let mp = MultiPolygon(vec![poly1, poly2]);
    /// // A larger single rectangle
    /// let combined = mp.unary_union();
    /// assert_eq!(&combined, merged);
    /// ```
    fn unary_union(self) -> MultiPolygon<Self::Scalar>;
}

// This function carries out a full post-order traversal of the tree, building up MultiPolygons from inside to outside.
// Though the operation is carried out via fold() over the tree iterator, there are two actual nested operations:
// "fold" operations on leaf nodes build up output MultiPolygons by adding Polygons to them via union and
// "reduce" operations on parent nodes combine these output MultiPolygons from leaf operations by recursion
fn bottom_up_fold_reduce<T, S, I, F, R>(tree: &RTree<T>, init: I, fold: F, reduce: R) -> S
where
    T: RTreeObject,
    RTreeNode<T>: Send + Sync,
    I: FnMut() -> S + Send + Sync,
    F: FnMut(S, &T) -> S + Send + Sync,
    R: FnMut(S, S) -> S + Send + Sync,
{
    // recursive algorithms can benefit from grouping those parameters which are constant over
    // the whole algorithm to reduce the overhead of the recursive calls
    struct Ops<I, F, R> {
        init: I,
        fold: F,
        reduce: R,
    }

    fn inner<T, S, I, F, R>(ops: &mut Ops<I, F, R>, parent: &ParentNode<T>) -> S
    where
        T: RTreeObject,
        I: FnMut() -> S,
        F: FnMut(S, &T) -> S,
        R: FnMut(S, S) -> S,
    {
        parent
            .children()
            .iter()
            .fold((ops.init)(), |accum, child| match child {
                RTreeNode::Leaf(value) => (ops.fold)(accum, value),
                RTreeNode::Parent(parent) => {
                    let value = inner(ops, parent);
                    (ops.reduce)(accum, value)
                }
            })
    }
    let mut ops = Ops { init, fold, reduce };
    inner(&mut ops, tree.root())
}

impl<T: BoolOpsNum> BooleanOps for Polygon<T> {
    type Scalar = T;

    fn rings(&self) -> impl Iterator<Item = &LineString<Self::Scalar>> {
        std::iter::once(self.exterior()).chain(self.interiors().iter())
    }
}

impl<T: BoolOpsNum> BooleanOps for MultiPolygon<T> {
    type Scalar = T;

    fn rings(&self) -> impl Iterator<Item = &LineString<Self::Scalar>> {
        self.0.iter().flat_map(|p| p.rings())
    }
}

// This struct and its RTReeObject impl allow construction of an R tree containing short-lived
// references to the original objects.
struct RTreeObjectRef<'a, T>(&'a T);

impl<'a, T> RTreeObject for RTreeObjectRef<'a, T>
where
    T: RTreeObject,
{
    type Envelope = T::Envelope;

    fn envelope(&self) -> Self::Envelope {
        self.0.envelope()
    }
}

impl<'a, T, Boppable, BoppableCollection> UnaryUnion for &'a BoppableCollection
where
    T: BoolOpsNum,
    Boppable: BooleanOps<Scalar = T> + RTreeObject + 'a + Sync,
    <Boppable as RTreeObject>::Envelope: Send + Sync,
    Self: IntoIterator<Item = &'a Boppable>,
{
    type Scalar = T;

    fn unary_union(self) -> MultiPolygon<Self::Scalar> {
        // these three functions drive the union operation
        let init = || MultiPolygon::<T>::new(vec![]);
        let fold = |mut accum: MultiPolygon<T>,
                    poly: &CachedEnvelope<RTreeObjectRef<'a, Boppable>>|
         -> MultiPolygon<T> {
            accum = accum.union(poly.0);
            accum
        };
        let reduce = |accum1: MultiPolygon<T>, accum2: MultiPolygon<T>| -> MultiPolygon<T> {
            accum1.union(&accum2)
        };
        let rtree = RTree::bulk_load(
            self.into_iter()
                .map(|p| CachedEnvelope::new(RTreeObjectRef(p)))
                .collect(),
        );
        bottom_up_fold_reduce(&rtree, init, fold, reduce)
    }
}

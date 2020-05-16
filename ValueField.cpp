#include "ValueField.h"
#include "Solver.h"

BaseValueField::BaseValueField(size_t str_size) {
	if (str_size > 0) resize_storage(str_size);
}

double VolumeBathymetryWrapper::at(index t) const {
	auto const& tp = b_.mesh().triang_points(t);
	return 1. / 3. * (b_.at_node(tp[0]) + b_.at_node(tp[1]) + b_.at_node(tp[2]));
};
constexpr index VolumeIndexer::Id(index i) {
	return i;
}
template struct ValueField<VolumeBathymetryWrapper, VolumeIndexer, index>;

double EdgeBathymetryWrapper::at(index edgeId, index fromId, index toId) const {
	auto const& ep = b_.mesh().edge_points(edgeId);
	return 0.5 * (b_.at_node(ep[0]) + b_.at_node(ep[1]));
};
constexpr index EdgeIndexer::Id(index edgeId, index fromId, index toId) {
	if (fromId == toId)
		throw SolverError("Equal ids of neighbours in edge indexer not supported");
	return 2 * edgeId + (index)(fromId < toId);
}
template struct ValueField<EdgeBathymetryWrapper, EdgeIndexer, index, index, index>;
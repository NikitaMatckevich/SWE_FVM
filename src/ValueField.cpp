#include <ValueField.h>
//#include <Solver.h>

BaseValueField::BaseValueField(size_t str_size) {
	if (str_size > 0) {
		resize_storage(str_size);
	}
} 

double VolumeBathymetryWrapper::at(idx t) const {
	auto const& tp = b_.mesh().triang_points(t);
	return 1. / 3. * (b_.at_node(tp[0]) + b_.at_node(tp[1]) + b_.at_node(tp[2]));
}
constexpr idx VolumeIndexer::Id(idx i) {
	return i;
}
template struct ValueField<VolumeBathymetryWrapper, VolumeIndexer, idx>;

double EdgeBathymetryWrapper::at(idx edgeId, idx fromId, idx toId) const {
	auto const& ep = b_.mesh().edge_points(edgeId);
	return 0.5 * (b_.at_node(ep[0]) + b_.at_node(ep[1]));
}
constexpr idx EdgeIndexer::Id(idx edgeId, idx fromId, idx toId) {
	if (fromId == toId)
		throw std::runtime_error("Equal ids of neighbours in edge indexer are not supported");
	return 2 * edgeId + (idx)(fromId < toId);
}
template struct ValueField<EdgeBathymetryWrapper, EdgeIndexer, idx, idx, idx>;

#include <ValueField.h>

BaseValueField::BaseValueField(size_t size) {
	if (size > 0)
		Resize(size);
}

double VolumeBathymetryWrapper::At(Idx t) const {
	const auto& tp = m_b.Mesh().TriangPoints(t);
	return (1. / 3.) * (m_b.AtNode(tp[0]) + m_b.AtNode(tp[1]) + m_b.AtNode(tp[2]));
}
template struct ValueField<VolumeBathymetryWrapper, VolumeIndexer, Idx>;

double EdgeBathymetryWrapper::At(Idx edgeId, Idx fromId, Idx toId) const {
	const auto& ep = m_b.Mesh().EdgePoints(edgeId);
	return 0.5 * (m_b.AtNode(ep[0]) + m_b.AtNode(ep[1]));
}
template struct ValueField<EdgeBathymetryWrapper, EdgeIndexer, Idx, Idx, Idx>;

#include <ValueField.h>

BaseValueField::BaseValueField(size_t size) {
	if (size > 0)
		Resize(size);
}

double VolumeDomainWrapper::At(Idx t) const {
	return m_b.T(t)[2];
}
template struct ValueField<VolumeDomainWrapper, VolumeIndexer, Idx>;

double EdgeDomainWrapper::At(Idx edgeId, Idx fromId, Idx toId) const {
	return m_b.E(edgeId)[2];
}
template struct ValueField<EdgeDomainWrapper, EdgeIndexer, Idx, Idx, Idx>;

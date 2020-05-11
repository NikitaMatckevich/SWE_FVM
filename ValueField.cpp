#include "ValueField.h"
#include "Solver.h"

namespace utils {

	BaseValueField::BaseValueField(const BaseValueField& other)
		: str_(other.str_), b_(other.b_) {};
	BaseValueField::BaseValueField(const Bathymetry& b, size_t str_sz)
		: b_(b)
	{
		str_.resize(Eigen::NoChange, str_sz);
	}
	const Bathymetry& BaseValueField::bathymetry() const { return b_; };

	constexpr index VolumeIndexer::Id(index i) {
		return i;
	}
	template struct ValueField<VolumeIndexer, index>;

	constexpr index EdgeIndexer::Id(index edgeId, index fromId, index toId) {
		if (fromId == toId)
			throw SolverError("Equal ids of neighbours in edge indexer not supported");
		return 2 * edgeId + (index)(fromId < toId);
	}
	template struct ValueField<EdgeIndexer, index, index, index>;

} // namespace utils
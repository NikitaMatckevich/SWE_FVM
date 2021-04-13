#pragma once
#include <ConsAssigner.h>
#include <iostream>
#include <cassert>

struct BaseValueField {

	explicit BaseValueField(size_t str_size = 0);
	
	size_t storage_size() const { return str_.cols(); }
	void resize_storage(size_t str_size) { str_.resize(Eigen::NoChange, str_size); }
  
 protected:
	Storage<3> str_;
};

template <class BathymetryWrapper, class Indexer, class ...Args>
struct ValueField : BaseValueField {

	const BathymetryWrapper b_;

	explicit ValueField(const Bathymetry& b, size_t str_size = 0)
		: BaseValueField(str_size), b_(b) {} 
	ValueField(const Bathymetry& b, const BaseValueField& other)
		: BaseValueField(other), b_(b) {}
	ValueField(const Bathymetry& b, BaseValueField&& other)
		: BaseValueField(std::move(other)), b_(b) {}

	auto     prim(Args... args) const { return str_.col(Indexer::Id(args...)); }
	Array<3> cons(Args... args) const { return { h(args...), hu(args...), hv(args...) }; }
	double w(Args... args) const { return str_(0, Indexer::Id(args...)); }
	double u(Args... args) const { return str_(1, Indexer::Id(args...)); }
	double v(Args... args) const { return str_(2, Indexer::Id(args...)); }
	auto vel (Args... args) const { return str_.col(Indexer::Id(args...)).tail(2); }
	double h (Args... args) const { return str_(0, Indexer::Id(args...)) - b_.at(args...); }
	double hu(Args... args) const { return h(args...) * u(args...); }
	double hv(Args... args) const { return h(args...) * v(args...); }

	PrimAssigner prim(Args... args) { return PrimAssigner(&str_, b_.at(args...), Indexer::Id(args...)); }
	ConsAssigner cons(Args... args) { return ConsAssigner(&str_, b_.at(args...), Indexer::Id(args...)); }
	double& u(Args... args) { return str_(1, Indexer::Id(args...)); }
	double& v(Args... args) { return str_(2, Indexer::Id(args...)); }
	auto  vel(Args... args) { return str_.col(Indexer::Id(args...)).tail(2); }
};

struct VolumeBathymetryWrapper : BaseBathymetryWrapper {
	using BaseBathymetryWrapper::BaseBathymetryWrapper;
	double at(idx t) const;
};
struct VolumeIndexer {
	static constexpr inline idx Id(idx t) {
		return t;
	}
};
using VolumeField = ValueField<VolumeBathymetryWrapper, VolumeIndexer, idx>;
extern template struct ValueField<VolumeBathymetryWrapper, VolumeIndexer, idx>;

struct EdgeBathymetryWrapper : BaseBathymetryWrapper {
	using BaseBathymetryWrapper::BaseBathymetryWrapper;
	double at(idx edgeId, idx fromId, idx toId) const;
} ;
struct EdgeIndexer {
	static constexpr inline idx Id(idx edgeId, idx fromId, idx toId) {
		assert(fromId != toId);
		return 2 * edgeId + (idx)(fromId < toId);
	}
};
using EdgeField = ValueField<EdgeBathymetryWrapper, EdgeIndexer, idx, idx, idx>;
extern template struct ValueField<EdgeBathymetryWrapper, EdgeIndexer, idx, idx, idx>;

#pragma once
#include "ConsAssigner.h"

struct BaseValueField {
	using Storage = Eigen::Array3Xd;
	using Array = Eigen::Array3d;
	explicit BaseValueField(size_t str_size = 0);
	size_t storage_size() const { return str_.cols(); }
	void resize_storage(size_t str_size) { str_.resize(Eigen::NoChange, str_size); }
protected:
	Storage str_;
};

template <class BathymetryWrapper, class Indexer, class ...Args>
struct ValueField : BaseValueField {
	const BathymetryWrapper b_;
	explicit ValueField(const Bathymetry& b, size_t str_size = 0)
		: BaseValueField(str_size), b_(b) {}

	auto  prim(Args... args) const { return str_.col(Indexer::Id(args...)); }
	Array cons(Args... args) const { return { h(args...), hu(args...), hv(args...) }; }
	double w(Args... args) const { return str_(0, Indexer::Id(args...)); }
	double u(Args... args) const { return str_(1, Indexer::Id(args...)); }
	double v(Args... args) const { return str_(2, Indexer::Id(args...)); }
	auto vel(Args... args) const { return str_.col(Indexer::Id(args...)).tail<2>(); }
	double h(Args... args) const { return str_(0, Indexer::Id(args...)) + b_.at(args...); }
	double hu(Args... args) const { return h(args...) * u(args...); }
	double hv(Args... args) const { return h(args...) * v(args...); }

	auto prim(Args... args) { return str_.col(Indexer::Id(args...)); }
	ConsAssigner cons(Args... args) { return ConsAssigner(&str_, b_.at(args...), Indexer::Id(args...)); };
	double& w(Args... args) { return str_(0, Indexer::Id(args...)); }
	double& u(Args... args) { return str_(1, Indexer::Id(args...)); }
	double& v(Args... args) { return str_(2, Indexer::Id(args...)); }
	auto vel(Args... args) { return str_.col(Indexer::Id(args...)).tail<2>(); }
};

struct VolumeBathymetryWrapper : BaseBathymetryWrapper {
	using BaseBathymetryWrapper::BaseBathymetryWrapper;
	double at(index t) const;
};
struct VolumeIndexer {
	static constexpr index Id(index t);
};
using VolumeField = ValueField<VolumeBathymetryWrapper, VolumeIndexer, index>;
extern template struct ValueField<VolumeBathymetryWrapper, VolumeIndexer, index>;

struct EdgeBathymetryWrapper : BaseBathymetryWrapper {
	using BaseBathymetryWrapper::BaseBathymetryWrapper;
	double at(index edgeId, index fromId, index toId) const;
};
struct EdgeIndexer {
	static constexpr index Id(index edgeId, index fromId, index toId);
};
using EdgeField = ValueField<EdgeBathymetryWrapper, EdgeIndexer, index, index, index>;
extern template struct ValueField<EdgeBathymetryWrapper, EdgeIndexer, index, index, index>;
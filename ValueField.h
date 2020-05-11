#pragma once
#include "ConsAssigner.h"

namespace utils {

	struct BaseValueField {
		using Storage = Eigen::Array3Xd;
		using Array = Eigen::Array3d;
		using ConstColXpr = Eigen::Block<const Storage, 3, 1, true>;
		BaseValueField() = default;
		BaseValueField(BaseValueField&& other) = default;
		BaseValueField(const BaseValueField& other);
		BaseValueField(const Bathymetry& b, size_t storage_size);
		const Bathymetry& bathymetry() const;
	protected:
		Storage str_;
	private:
		const Bathymetry& b_;
	};

	template <class Indexer, class ...Args>
	struct ValueField : BaseValueField {
		using Base = BaseValueField;
		using Array = Base::Array;
		using Base::BaseValueField;

		auto prim(Args... args) const { return str_.col(Indexer::Id(args...)); }
		Array cons(Args... args) const { return { h(args...), hu(args...), hv(args...) }; }
		double w(Args... args) const { return str_(0, Indexer::Id(args...)); }
		double u(Args... args) const { return str_(1, Indexer::Id(args...)); }
		double v(Args... args) const { return str_(2, Indexer::Id(args...)); }
		auto vel(Args... args) const { return str_.col(Indexer::Id(args...)).tail<2>(); }
		double h(Args... args) const {
			auto i = Indexer::Id(args...);
			return str_(0, i) + bathymetry().t(i);
		}
		double hu(Args... args) const { return h(args...) * u(args...); }
		double hv(Args... args) const { return h(args...) * v(args...); }

		auto prim(Args... args) { return str_.col(Indexer::Id(args...)); }
		ConsAssigner cons(Args... args) { return ConsAssigner(&str_, bathymetry(), Indexer::Id(args...)); };
		double& w(Args... args) { return str_(0, Indexer::Id(args...)); }
		double& u(Args... args) { return str_(1, Indexer::Id(args...)); }
		double& v(Args... args) { return str_(2, Indexer::Id(args...)); }
		auto vel(Args... args) { return str_.col(Indexer::Id(args...)).tail<2>(); }
	};

	struct VolumeIndexer {
		static constexpr index Id(index i);
	};
	using VolumeField = ValueField<VolumeIndexer, index>;
	extern template struct ValueField<VolumeIndexer, index>;

	struct EdgeIndexer {
		static constexpr index Id(index edgeId, index fromId, index toId);
	};
	using EdgeField = ValueField<EdgeIndexer, index, index, index>;
	extern template struct ValueField<EdgeIndexer, index, index, index>;

} // namespace utils
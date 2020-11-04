#pragma once
#include <fluxions_gte.hpp>

namespace Sf {
	namespace Fx = Fluxions;

	class IPBSky {
	public:
		using Real = float;
		IPBSky() {}
		virtual ~IPBSky() {}

		virtual void init(Real turbidity, Fx::Color4f albedo_, Real elevation, Real azimuth) = 0;
		virtual void kill() = 0;
		virtual void resetStatisticSamples() = 0;
		virtual void addStatisticSample(Real amount) const = 0;

		virtual void computeThetaGamma(Real inclination, Real azimuth, Real* outTheta, Real* outGamma) const = 0;
		virtual void computeThetaGamma(Real x, Real y, Real z, Real* outTheta, Real* outGamma) const = 0;
		virtual void computeSunRadiance4(Real theta, Real gamma, Fx::Color4f& output) const = 0;

		virtual Fx::Color4f getSunDiskRadiance() const = 0;
		virtual Fx::Color4f getGroundRadiance() const = 0;

		virtual Real getAverageRadiance() const = 0;
		virtual Real getMinValue() const = 0;
		virtual Real getMaxValue() const = 0;
	};
}

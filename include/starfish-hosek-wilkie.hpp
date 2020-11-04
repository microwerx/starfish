#pragma once
#include <ArHosekSkyModel.h>
#include <fluxions_gte.hpp>
#include <starfish-ipbsky.hpp>

namespace Sf {
	namespace Fx = Fluxions;

	class HosekWilkiePBSky : public IPBSky {
	public:
		using Real = float;

		HosekWilkiePBSky();
		HosekWilkiePBSky(Real turbidity, Fx::Color4f albedo, Real elevation, Real azimuth);
		~HosekWilkiePBSky() override;

		void init(Real turbidity, Fx::Color4f albedo_, Real elevation, Real azimuth) override;
		void kill() override;
		void resetStatisticSamples() override;
		void addStatisticSample(Real amount) const override;

		void computeThetaGamma(Real inclination, Real azimuth, Real* outTheta, Real* outGamma) const override;
		void computeThetaGamma(Real x, Real y, Real z, Real* outTheta, Real* outGamma) const override;
		void computeSunRadiance4(Real theta, Real gamma, Fx::Color4f& output) const override;

		Fx::Color4f getSunDiskRadiance() const override;
		Fx::Color4f getGroundRadiance() const override;

		Real getAverageRadiance() const override {
			return totalValue / Real(nSamples);
		}

		Real getMinValue() const override {
			return minValue;
		}

		Real getMaxValue() const override {
			return maxValue;
		}

	private:
		Real sunAltitude{ 0 };
		Real sunInclination{ 0 };
		Real sunAzimuth{ 0 };
		Real sun[3]{ 0,0,0 };
		Real sunTheta{ 0 };
		Real sunGamma{ 0 };
		Real sunTurbidity{ 1 };
		Real sunElevation{ 0 };
		Fx::Color4f albedo;

		Real compute(Real theta, Real gamma, int index);
		void computeSunRadiance(Real theta, Real gamma, Fx::Color4f& output) const;
		Real computeSunRadiance2(Real theta, Real gamma, int index);
		// 11 band Solar Radiance
		void computeSunRadiance3(Real theta, Real gamma, Fx::Color4f& output);
		// 3 band RGB Radiance
		void computeSunRadiance4_NoStatistics(Real theta, Real gamma, Fx::Color4f& output) const;

	private:
		ArHosekSkyModelState* rgbRadianceState[3];
		ArHosekSkyModelState* sunRadianceState; // 320nm to 720nm in steps of 40nm for wavelength

		mutable volatile Real minValue{ 0 };
		mutable volatile Real maxValue{ 0 };
		mutable volatile Real totalValue{ 0 };
		mutable volatile int nSamples{ 0 };
	};
} // namespace Sf

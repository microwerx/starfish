#include "starfish_pch.hpp"
#include <starfish-hosek-wilkie.hpp>
#include <hatchetfish.hpp>

namespace Sf {
	using Real = float;
	const int NUM_WAVELENGTHS = 11;


	Real wavelengths[NUM_WAVELENGTHS] = {
		320.0f,
		360.0f,
		400.0f,
		440.0f,
		480.0f,
		520.0f,
		560.0f,
		600.0f,
		640.0f,
		680.0f,
		720.0f };


	// weights to integrate the tristimulus values about a
	Real tristimulus[NUM_WAVELENGTHS][3] = {
		{0.0f, 0.0f, 0.0f},
		{0.0f, 0.0f, 0.0f},
		{0.0174f, 0.0018f, 0.0881f},
		{0.1056f, 0.0173f, 0.585f},
		{0.037f, 0.0656f, 0.2973f},
		{0.0359f, 0.2406f, 0.0287f},
		{0.2179f, 0.3455f, 0.0009f},
		{0.3822f, 0.2433f, 0.0f},
		{0.1805f, 0.0768f, 0.0f},
		{0.0224f, 0.0087f, 0.0f},
		{0.0013f, 0.0005f, 0.0f} };


	HosekWilkiePBSky::HosekWilkiePBSky() {
		for (int i = 0; i < 3; i++)
			rgbRadianceState[i] = nullptr;
		sunRadianceState = nullptr;
		//for (int i = 0; i < NUM_WAVELENGTHS; i++) sunRadianceState[i] = nullptr;
	}


	HosekWilkiePBSky::HosekWilkiePBSky(Real turbidity, Fx::Color4f albedo, Real elevation, Real azimuth) {
		for (int i = 0; i < 3; i++)
			rgbRadianceState[i] = nullptr;
		sunRadianceState = nullptr;
		//for (int i = 0; i < 11; i++) sunRadianceState[i] = nullptr;
		init(turbidity, albedo, elevation, azimuth);
	}


	HosekWilkiePBSky::~HosekWilkiePBSky() {
		// Delete();
	}


	void HosekWilkiePBSky::init(Real turbidity, Fx::Color4f albedo_, Real elevation, Real azimuth) {
		// offset by -90 degrees because we are measuring from NORTH which is positive Y
		azimuth = 90.0f - azimuth;
		elevation *= (float)FX_DEGREES_TO_RADIANS;
		azimuth *= (float)FX_DEGREES_TO_RADIANS;

		if (this->albedo.r != albedo_.r || this->albedo.g != albedo_.g || this->albedo.b != albedo_.b || this->sunTurbidity != turbidity || this->sunElevation != elevation) {
			kill();
		}

		if (!rgbRadianceState[0])
			rgbRadianceState[0] = arhosek_rgb_skymodelstate_alloc_init(turbidity, albedo_.r, elevation);
		if (!rgbRadianceState[1])
			rgbRadianceState[1] = arhosek_rgb_skymodelstate_alloc_init(turbidity, albedo_.g, elevation);
		if (!rgbRadianceState[2])
			rgbRadianceState[2] = arhosek_rgb_skymodelstate_alloc_init(turbidity, albedo_.b, elevation);

		sunRadianceState = arhosekskymodelstate_alloc_init(elevation, turbidity, albedo_.r);
		//for (int i = 0; i < NUM_WAVELENGTHS; i++)
		//{
		//	sunRadianceState[i] = arhosekskymodelstate_alloc_init(elevation, turbidity, albedo_.r);
		//}

		this->albedo = albedo_;
		sunTurbidity = turbidity;
		sunElevation = elevation;

		minValue = FLT_MAX;
		maxValue = -FLT_MAX;
		totalValue = 0.0;
		nSamples = 0;

		sunAltitude = elevation;
		sunAzimuth = azimuth;
		sunInclination = (float)FX_PIOVERTWO - sunAltitude;
		sun[0] = sin(sunInclination) * cos(sunAzimuth);
		sun[1] = sin(sunInclination) * sin(sunAzimuth);
		sun[2] = cos(sunInclination);
		computeThetaGamma(sun[0], sun[1], sun[2], &sunTheta, &sunGamma);
	}


	void HosekWilkiePBSky::kill() {
		try {
			for (int i = 0; i < 3; i++) {
				if (rgbRadianceState[i])
					arhosekskymodelstate_free(rgbRadianceState[i]);
				rgbRadianceState[i] = nullptr;
			}

			if (sunRadianceState)
				arhosekskymodelstate_free(sunRadianceState);
			sunRadianceState = nullptr;
			//for (int i = 0; i < NUM_WAVELENGTHS; i++)
			//{
			//	if (sunRadianceState[i])
			//		arhosekskymodelstate_free(sunRadianceState[i]);
			//	sunRadianceState[i] = nullptr;
			//}
		}
		catch (...) {
			HFLOGERROR("unknown error freeing memory.");
		}
	}


	void HosekWilkiePBSky::resetStatisticSamples() {
		minValue = 1e10;
		maxValue = -1e10;
		totalValue = 0.0;
		nSamples = 0;
	}


	void HosekWilkiePBSky::addStatisticSample(Real amount) const {
		if (minValue > amount) {
			minValue = amount;
		}
		if (maxValue < amount) {
			maxValue = amount;
		}
		totalValue += amount;
		nSamples++;
	}


	void HosekWilkiePBSky::computeThetaGamma(Real inclination, Real azimuth, Real* outTheta, Real* outGamma) const {
		Real v[3];
		v[0] = sin(inclination) * cos(azimuth);
		v[1] = sin(inclination) * sin(azimuth);
		v[2] = cos(inclination);
		// Gamma is angle between sun vector and sky element vector which is the dot product
		Real cosine = v[0] * sun[0] + v[1] * sun[1] + v[2] * sun[2];
		*outGamma = acos(cosine);
		*outTheta = inclination;
	}


	void HosekWilkiePBSky::computeThetaGamma(Real x, Real y, Real z, Real* outTheta, Real* outGamma) const {
		Real length = sqrtf(x * x + y * y + z * z);
		Real v[3];
		v[0] = x / length;
		v[1] = y / length;
		v[2] = z / length;
		Real cosine = std::clamp(v[0] * sun[0] + v[1] * sun[1] + v[2] * sun[2], -1.0f, 1.0f);
		Real gamma = acos(cosine);
		Real theta = acos(v[2]);
		if (std::isfinite(gamma)) {
			*outGamma = gamma;
		}
		else {
			*outGamma = 0.0;
		}
		if (std::isfinite(theta)) {
			*outTheta = theta;
		}
		else {
			*outTheta = 0.0;
		}
	}


	Real HosekWilkiePBSky::compute(Real theta, Real gamma, int index) {
		if (index < 0 || index >= 3)
			return 0.0;

		if (theta < 0.0 || theta >= FX_PIOVERTWO || gamma < 0.0)
			return 0.0;

		Real amount = (float)arhosek_tristim_skymodel_radiance(rgbRadianceState[index], theta, gamma, index);

		addStatisticSample(amount);

		return amount;
	}


	void HosekWilkiePBSky::computeSunRadiance(Real theta, Real gamma, Fx::Color4f& output) const {
		Fx::Color4f xyz;

		for (int i = 0; i < NUM_WAVELENGTHS; i++) {
			Real amount = theta > FX_PIOVERTWO ? 0.0f : (float)arhosekskymodel_solar_radiance(sunRadianceState, theta, gamma, wavelengths[i]);
			xyz.r += tristimulus[i][0] * amount;
			xyz.g += tristimulus[i][1] * amount;
			xyz.b += tristimulus[i][2] * amount;
		}

		// convert XYZ to sRGB
		Real m[3][3] = {
			{3.2406f, -1.5372f, -0.4986f},
			{-0.9689f, 1.8758f, 0.0415f},
			{0.0557f, -0.2040f, 1.0570f} };

		Fx::Color4f rgb_linear(
			xyz.r * m[0][0] + xyz.g * m[0][1] + xyz.b * m[0][2],
			xyz.r * m[1][0] + xyz.g * m[1][1] + xyz.b * m[1][2],
			xyz.r * m[2][0] + xyz.g * m[2][1] + xyz.b * m[2][2],
			0.0f);

		Fx::Color4f sRGB(
			rgb_linear.r <= 0.0031308f ? rgb_linear.r * 12.92f : (1 + 0.055f) * powf(rgb_linear.r, 1.0f / 2.4f) - 0.055f,
			rgb_linear.g <= 0.0031308f ? rgb_linear.g * 12.92f : (1 + 0.055f) * powf(rgb_linear.g, 1.0f / 2.4f) - 0.055f,
			rgb_linear.b <= 0.0031308f ? rgb_linear.b * 12.92f : (1 + 0.055f) * powf(rgb_linear.b, 1.0f / 2.4f) - 0.055f,
			0.0f);

		output = sRGB;
	}


	Real HosekWilkiePBSky::computeSunRadiance2(Real theta, Real gamma, int index) {
		//return Compute(theta, gamma, index);

		if (theta < 0.0 || theta >= FX_PIOVERTWO || gamma < 0.0)
			return 0.0;
		Real amount;

		Fx::Color4f xyz;

		for (int i = 0; i < NUM_WAVELENGTHS; i++) {
			//Real amount = arhosekskymodel_solar_radiance(sunRadianceState[i], theta, gamma, wavelengths[i]);
			amount = (float)arhosekskymodel_inscattered_radiance(sunRadianceState, theta, gamma, wavelengths[i]);
			if (theta < 1.0f)
				amount += (float)arhosekskymodel_sun_direct_radiance(sunRadianceState, theta, gamma, wavelengths[i]);
			//amount = (float)arhosekskymodel_solar_radiance(sunRadianceState[i], theta, gamma, wavelengths[i]);

			if (std::isfinite(amount))
				addStatisticSample(amount);
			//else std::cout << "BLAH!\n";

			xyz.r += tristimulus[i][0] * amount;
			xyz.g += tristimulus[i][1] * amount;
			xyz.b += tristimulus[i][2] * amount;
		}
		if (index == 0)
			return xyz.r;
		if (index == 1)
			return xyz.g;
		if (index == 2)
			return xyz.b;

		// convert XYZ to sRGB
		Real m[3][3] = {
			{3.2406f, -1.5372f, -0.4986f},
			{-0.9689f, 1.8758f, 0.0415f},
			{0.0557f, -0.2040f, 1.0570f} };

		Fx::Color4f rgb_linear(
			xyz.r * m[0][0] + xyz.g * m[0][1] + xyz.b * m[0][2],
			xyz.r * m[1][0] + xyz.g * m[1][1] + xyz.b * m[1][2],
			xyz.r * m[2][0] + xyz.g * m[2][1] + xyz.b * m[2][2],
			0.0f);

		Fx::Color4f sRGB(
			rgb_linear.r <= 0.0031308f ? rgb_linear.r * 12.92f : (1 + 0.055f) * powf(rgb_linear.r, 1.0f / 2.4f) - 0.055f,
			rgb_linear.g <= 0.0031308f ? rgb_linear.g * 12.92f : (1 + 0.055f) * powf(rgb_linear.g, 1.0f / 2.4f) - 0.055f,
			rgb_linear.b <= 0.0031308f ? rgb_linear.b * 12.92f : (1 + 0.055f) * powf(rgb_linear.b, 1.0f / 2.4f) - 0.055f,
			0.0f);

		if (index == 0)
			return sRGB.r;
		if (index == 1)
			return sRGB.g;
		if (index == 2)
			return sRGB.b;
		return 0.0f;
	}


	void HosekWilkiePBSky::computeSunRadiance3(Real theta, Real gamma, Fx::Color4f& output) {
		//return Compute(theta, gamma, index);

		if (theta < 0.0 || theta >= FX_PIOVERTWO || gamma < 0.0)
			return;
		Real amount;

		Fx::Color4f xyz;

		for (int i = 0; i < NUM_WAVELENGTHS; i++) {
			//Real amount = arhosekskymodel_solar_radiance(sunRadianceState[i], theta, gamma, wavelengths[i]);
			amount = (float)arhosekskymodel_inscattered_radiance(sunRadianceState, theta, gamma, wavelengths[i]);

			if (std::isfinite(amount))
				addStatisticSample(amount);

			if (gamma < 0.01f)
				amount += (float)arhosekskymodel_sun_direct_radiance(sunRadianceState, theta, gamma, wavelengths[i]);
			//amount = (float)arhosekskymodel_solar_radiance(sunRadianceState[i], theta, gamma, wavelengths[i]);

			//else std::cout << "BLAH!\n";

			xyz.r += tristimulus[i][0] * amount;
			xyz.g += tristimulus[i][1] * amount;
			xyz.b += tristimulus[i][2] * amount;
		}
		output.r = xyz.r;
		output.g = xyz.g;
		output.b = xyz.b;
	}


	void HosekWilkiePBSky::computeSunRadiance4(Real theta, Real gamma, Fx::Color4f& output) const {
		Fx::Color4f xyz;

		// Return 0 if we are outside of the domain of this function.
		if (theta < 0.0 || theta >= FX_PIOVERTWO || gamma < 0.0) {
			output.r = 0.5f;
			output.g = 0.5f;
			output.b = 0.5f;
			output.a = 1.0f;
			return;
		}

		output.r = (float)arhosek_tristim_skymodel_radiance(rgbRadianceState[0], theta, gamma, 0);
		output.g = (float)arhosek_tristim_skymodel_radiance(rgbRadianceState[1], theta, gamma, 1);
		output.b = (float)arhosek_tristim_skymodel_radiance(rgbRadianceState[2], theta, gamma, 2);
		output.a = 1.0f;

		if (std::isfinite(output.r))
			addStatisticSample(output.r);
		if (std::isfinite(output.g))
			addStatisticSample(output.g);
		if (std::isfinite(output.b))
			addStatisticSample(output.b);
	}


	void HosekWilkiePBSky::computeSunRadiance4_NoStatistics(Real theta, Real gamma, Fx::Color4f& output) const {
		Fx::Color4f xyz;

		// Return 0 if we are outside of the domain of this function.
		if (theta < 0.0 || theta >= FX_PIOVERTWO || gamma < 0.0) {
			output.r = 0.0f;
			output.g = 0.0f;
			output.b = 0.0f;
			output.a = 1.0f;
			return;
		}

		output.r = (float)arhosek_tristim_skymodel_radiance(rgbRadianceState[0], theta, gamma, 0);
		output.g = (float)arhosek_tristim_skymodel_radiance(rgbRadianceState[1], theta, gamma, 1);
		output.b = (float)arhosek_tristim_skymodel_radiance(rgbRadianceState[2], theta, gamma, 2);
		output.a = output.ToVector3().length();
	}


	Fx::Color4f HosekWilkiePBSky::getSunDiskRadiance() const {
		Fx::Color4f output;
		//ComputeSunRadiance4_NoStatistics(sunTheta, sunGamma, output);
		computeSunRadiance(sunTheta, sunGamma, output);

		if (!std::isnormal(output.r) || output.r < 0.0f)
			output.r = 0.0f;
		if (!std::isnormal(output.g) || output.g < 0.0f)
			output.g = 0.0f;
		if (!std::isnormal(output.b) || output.b < 0.0f)
			output.b = 0.0f;
		if (!std::isnormal(output.a) || output.a < 0.0f)
			output.a = 0.0f;

		return Fx::Color4f(output);
	}


	Fx::Color4f HosekWilkiePBSky::getGroundRadiance() const {
		Fx::Color4f sunRadiance = getSunDiskRadiance();
		float f_r = (float)(1.0 / FX_PI * std::max(0.0f, sun[2]));

		return Fx::Color4f(
			sunRadiance.r * albedo.r * f_r,
			sunRadiance.g * albedo.g * f_r,
			sunRadiance.b * albedo.b * f_r,
			1.0);
	}
}

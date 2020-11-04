#include "starfish_pch.hpp"
#include <hatchetfish.hpp>
#include <starfish-pbsky.hpp>
#include <fluxions_base.hpp>
#include <fluxions_image_loader.hpp>


namespace Sf {
	namespace Fx = Fluxions;
	using Fx::Vector3f;
	using Fx::Color3f;
	using Fx::Color4f;
	using Fx::Image3f;



	// Determines which face should be selected for rendering/getting from
	int classifyCubeFaceFromVector(float x, float y, float z) {
		// ma is absolute value
		float ax = fabs(x);
		float ay = fabs(y);
		float az = fabs(z);
		// signs of ma
		int sx = x > 0 ? 1 : 0;
		int sy = y > 0 ? 1 : 0;
		int sz = z > 0 ? 1 : 0;

		// GL_TEXTURE_CUBE_MAP_POSITIVE_X
		if (sx && ax >= ay && ax >= az)
			return 0;
		// GL_TEXTURE_CUBE_MAP_NEGATIVE_X
		if (!sx && ax >= ay && ax >= az)
			return 1;
		// GL_TEXTURE_CUBE_MAP_POSITIVE_Y
		if (sy && ay >= ax && ay >= az)
			return 2;
		// GL_TEXTURE_CUBE_MAP_NEGATIVE_Y
		if (!sy && ay >= ax && ay >= az)
			return 3;
		// GL_TEXTURE_CUBE_MAP_POSITIVE_Z
		if (sz && az >= ax && az >= ay)
			return 4;
		// GL_TEXTURE_CUBE_MAP_NEGATIVE_Z
		if (!sz && az >= ax && az >= ay)
			return 5;

		return -1;
	}

	void makeST(float x, float y, float z, float* s, float* t, int* whichFace) {
		// ma is absolute value
		float ax = fabs(x);
		float ay = fabs(y);
		float az = fabs(z);
		// signs of ma
		int sx = x > 0 ? 1 : 0;
		int sy = y > 0 ? 1 : 0;
		int sz = z > 0 ? 1 : 0;

		// 0 | sc = -z | tc = -y | ma = +X | dir: +rx
		// 1 | sc = +z | tc = -y | ma = +X | dir: -rx
		// 2 | sc = +X | tc = +z | ma = +y | dir: +ry
		// 3 | sc = +X | tc = -z | ma = +y | dir: -ry
		// 4 | sc = +X | tc = -y | ma = +z | dir: +rz
		// 5 | sc = -X | tc = -y | ma = +z | dir: -rz
		float ma = 0.0f;
		float sc = 0.0f;
		float tc = 0.0f;
		// GL_TEXTURE_CUBE_MAP_POSITIVE_X
		if (sx && ax >= ay && ax >= az) {
			ma = ax;
			sc = z;
			tc = y;
			*whichFace = 0;
		}
		// GL_TEXTURE_CUBE_MAP_NEGATIVE_X
		if (!sx && ax >= ay && ax >= az) {
			ma = ax;
			sc = z;
			tc = y;
			*whichFace = 1;
		}
		// GL_TEXTURE_CUBE_MAP_POSITIVE_Y
		if (sy && ay >= ax && ay >= az) {
			ma = ay;
			sc = x;
			tc = z;
			*whichFace = 2;
		}
		// GL_TEXTURE_CUBE_MAP_NEGATIVE_Y
		if (!sy && ay >= ax && ay >= az) {
			ma = ay;
			sc = x;
			tc = z;
			*whichFace = 3;
		}
		// GL_TEXTURE_CUBE_MAP_POSITIVE_Z
		if (sz && az >= ax && az >= ay) {
			ma = az;
			sc = x;
			tc = y;
			*whichFace = 4;
		}
		// GL_TEXTURE_CUBE_MAP_NEGATIVE_Z
		if (!sz && az >= ax && az >= ay) {
			ma = az;
			sc = x;
			tc = y;
			*whichFace = 5;
		}
		// s = 0.5 * (sc / ma + 1)
		// t = 0.5 * (tc / ma + 1)
		*s = 0.5f * (sc / ma + 1.0f);
		*t = 0.5f * (tc / ma + 1.0f);
	}

	Vector3f makeCubeVector(int face, float s, float t) {
		// 0 GL_TEXTURE_CUBE_MAP_POSITIVE_X
		// 1 GL_TEXTURE_CUBE_MAP_NEGATIVE_X
		// 2 GL_TEXTURE_CUBE_MAP_POSITIVE_Y
		// 3 GL_TEXTURE_CUBE_MAP_NEGATIVE_Y
		// 4 GL_TEXTURE_CUBE_MAP_POSITIVE_Z
		// 5 GL_TEXTURE_CUBE_MAP_NEGATIVE_Z
		//
		// ma is the major axis (the largest of the absolute values of X, y, z)
		// sc is the s coordinate in 3D
		// tc is the t coordinate in 3D
		// 0 | sc = -z | tc = -y | ma = +X | dir: +rx
		// 1 | sc = +z | tc = -y | ma = +X | dir: -rx
		// 2 | sc = +X | tc = +z | ma = +y | dir: +ry
		// 3 | sc = +X | tc = -z | ma = +y | dir: -ry
		// 4 | sc = +X | tc = -y | ma = +z | dir: +rz
		// 5 | sc = -X | tc = -y | ma = +z | dir: -rz
		//
		// s = 0.5 * (sc / ma + 1)
		// t = 0.5 * (tc / ma + 1)
		//
		// sc = ma (2s - 1)
		// tc = ma (2t - 1)
		//
		float x = 0.0f;
		float y = 0.0f;
		float z = 0.0f;
		float sc = 2.0f * s - 1;
		float tc = 2.0f * t - 1;

		switch (face) {
		case 0:
			x = 1.0;
			y = tc;
			z = -sc;
			break;
		case 1:
			x = -1.0;
			y = tc;
			z = sc;
			break;
		case 2:
			x = sc;
			y = 1.0;
			z = -tc;
			break;
		case 3:
			x = sc;
			y = -1.0;
			z = tc;
			break;
		case 4:
			x = sc;
			y = tc;
			z = 1.0;
			break;
		case 5:
			x = -sc;
			y = tc;
			z = -1.0;
			break;
		}

		return Vector3f(x, y, z);

		//// Test this by going the other way...
		//float scheck, tcheck;
		//int whichFace;
		//makeST(x, y, z, &scheck, &tcheck, &whichFace);

		//// check if face maps out. =^)
		////whichFace = classifyCubeFaceFromVector(X, y, z);

		//cout << std::setprecision(2) << left;
		//cout << "S: " << std::setw(5) << s << " T: " << std::setw(5) << t << " ";
		//cout << "S: " << std::setw(5) << scheck << " T: " << std::setw(5) << tcheck << " ";
		//cout << "Face: " << face << " ";
		//cout << "<--> " << whichFace << " ";
		//cout << showpos << left << showpoint;
		//cout << "(" << std::setw(7) << x << ", " << std::setw(7) << y << ", " << z << ") ";
		//cout << noshowpos;
		//cout << std::endl;

		//return Vector3f(x, y, z);
	}

	Vector3f makeCubeVector2(int face, float s, float t) {
		// 0 GL_TEXTURE_CUBE_MAP_POSITIVE_X
		// 1 GL_TEXTURE_CUBE_MAP_NEGATIVE_X
		// 2 GL_TEXTURE_CUBE_MAP_POSITIVE_Y
		// 3 GL_TEXTURE_CUBE_MAP_NEGATIVE_Y
		// 4 GL_TEXTURE_CUBE_MAP_POSITIVE_Z
		// 5 GL_TEXTURE_CUBE_MAP_NEGATIVE_Z
		//
		// ma is the major axis (the largest of the absolute values of X, y, z)
		// sc is the s coordinate in 3D
		// tc is the t coordinate in 3D
		// 0 | sc = -z | tc = -y | ma = +X | dir: +rx
		// 1 | sc = +z | tc = -y | ma = +X | dir: -rx
		// 2 | sc = +X | tc = +z | ma = +y | dir: +ry
		// 3 | sc = +X | tc = -z | ma = +y | dir: -ry
		// 4 | sc = +X | tc = -y | ma = +z | dir: +rz
		// 5 | sc = -X | tc = -y | ma = +z | dir: -rz
		//
		// s = 0.5 * (sc / ma + 1)
		// t = 0.5 * (tc / ma + 1)
		//
		// sc = ma (2s - 1)
		// tc = ma (2t - 1)
		//
		float x = 0.0f;
		float y = 0.0f;
		float z = 0.0f;
		float sc = 2.0f * s - 1;
		float tc = 2.0f * t - 1;

		switch (face) {
		case 0:
			x = 1.0;
			y = tc;
			z = -sc;
			break;
		case 1:
			x = -1.0;
			y = tc;
			z = sc;
			break;
		case 2:
			x = sc;
			y = 1.0;
			z = -tc;
			break;
		case 3:
			x = sc;
			y = -1.0;
			z = tc;
			break;
		case 4:
			x = sc;
			y = tc;
			z = 1.0;
			break;
		case 5:
			x = -sc;
			y = tc;
			z = -1.0;
			break;
		}

		return Vector3f(x, y, z);

		// TODO: Move this to a separate test

		//if ((rand() % 400) != 1)
		//	return Vector3f(x, y, z);

		//// Test this by going the other way...
		//float scheck, tcheck;
		//int whichFace;
		//makeST(x, y, z, &scheck, &tcheck, &whichFace);

		//// check if face maps out. =^)
		//whichFace = classifyCubeFaceFromVector(x, y, z);

		//std::cout << std::setprecision(2) << std::left;
		//std::cout << "S: " << std::setw(5) << s << " T: " << std::setw(5) << t << " ";
		//std::cout << "S: " << std::setw(5) << scheck << " T: " << std::setw(5) << tcheck << " ";
		//std::cout << "Face: " << face << " ";
		//if (whichFace >= 0)
		//	std::cout << "<--> " << whichFace << " ";
		//else
		//	std::cout << "<--> X ";
		//std::cout << std::showpos << std::left << std::showpoint;
		//std::cout << "(" << std::setw(7) << x << ", " << std::setw(7) << y << ", " << z << ") ";
		//std::cout << std::noshowpos;
		//std::cout << std::endl;

		//return Vector3f(x, y, z);
	}

	/////////////////////////////////////////////////////////////////////
	// P H Y S I C A L L Y   B A S E D   S K Y //////////////////////////
	/////////////////////////////////////////////////////////////////////

	PhysicallyBasedSky::PhysicallyBasedSky() {
		hwSunPbsky = std::make_shared<HosekWilkiePBSky>();
		hwMoonPbsky = std::make_shared<HosekWilkiePBSky>();
	}

	PhysicallyBasedSky::~PhysicallyBasedSky() {}

	void PhysicallyBasedSky::setLocation(float latitude, float longitude) {
		astroCalc.SetLocation(latitude, longitude);
		computeAstroFromLocale();
	}

	time_t PhysicallyBasedSky::getTime() const {
		return astroCalc.GetTime();
	}

	void PhysicallyBasedSky::setTime(time_t t, float fractSeconds) {
		astroCalc.SetTime(t, fractSeconds);
		computeAstroFromLocale();
	}

	void PhysicallyBasedSky::setLocalDate(int day, int month, int year, bool isdst, int timeOffset) {
		astroCalc.SetDate(day, month, year, isdst, timeOffset);
		computeAstroFromLocale();
	}

	void PhysicallyBasedSky::setCivilDateTime(const PA::CivilDateTime& dtg) {
		astroCalc.SetDateTime(dtg.day, dtg.month, dtg.year, dtg.isdst, dtg.timeZoneOffset, dtg.hh, dtg.mm, dtg.ss, dtg.ss_frac);
		computeAstroFromLocale();
	}

	void PhysicallyBasedSky::setLocalTime(int hh, int mm, int ss, float ss_frac) {
		astroCalc.SetTime(hh, mm, ss, ss_frac);
		computeAstroFromLocale();
	}

	void PhysicallyBasedSky::setTurbidity(float T) {
		turbidity = T;
	}

	float PhysicallyBasedSky::getTurbidity() const {
		return turbidity;
	}

	void PhysicallyBasedSky::setSunPosition(double azimuth, double altitude) {
		sunPosition_.A = azimuth;
		sunPosition_.a = altitude;
	}

	void PhysicallyBasedSky::setSunPosition(double sunLong) {
		EclipticCoord sunCoord(sunLong, 0.0);
		sunRADec_ = astroCalc.ecliptic_to_equatorial(sunCoord);
		sunPosition_ = astroCalc.ecliptic_to_horizon(sunCoord);
		Vector v = sunPosition_.toOpenGLVector();
		sunVector_ = { (float)v.x, (float)v.y, (float)v.z };
	}

	void PhysicallyBasedSky::setMoonPosition(double RA, double dec) {
		setMoonPosition({ RA, dec });
	}

	void PhysicallyBasedSky::setMoonPosition(EquatorialCoord moonRADec) {
		moonRADec_ = moonRADec;
		HorizonCoord p = astroCalc.equatorial_to_horizon(moonRADec);
		Vector v = p.toOpenGLVector();
		moonVector_ = { (float)v.x, (float)v.y, (float)v.z };
	}

	float PhysicallyBasedSky::getAverageRadiance() const {
		return hwSunPbsky->getAverageRadiance();
	}

	void PhysicallyBasedSky::computeAstroFromLocale() {
		_computeSunFromLocale();
		_computeMoonFromLocale();
	}

	void PhysicallyBasedSky::_computeSunFromLocale() {
		EclipticCoord sunEcl = astroCalc.getSun();
		double sunLong = Fx::wrap<double>(sunEcl.lambda, 360.0);
		setSunPosition(sunLong);
	}

	void PhysicallyBasedSky::_computeMoonFromLocale() {
		setMoonPosition(astroCalc.getMoon());
	}

	void PhysicallyBasedSky::computeCubeMap(int resolution, bool normalize, float sampleScale, bool flipY) {
		_prepareForCompute();
		generatedSunCubeMap.resize(resolution, resolution, 6);

		float sampleRadius = nSamples > 1 ? 0.5f / (float)resolution : 0.0f;

		// Loop through each texel and compute (s, t) coordinates such that each
		// ray goes the center of each pixel.
		for (int face = 0; face < 6; face++) {
			for (int is = 0; is < resolution; is++) {
				for (int it = 0; it < resolution; it++) {
					float s = (is + 0.5f) / (float)resolution;
					float t = (it + 0.5f) / (float)resolution;
					Color4f color;

					for (int curSample = 0; curSample < nSamples; curSample++) {
						float sp = s;
						float tp = t;
						if (curSample >= 1) {
							sp += Fx::randomSampler(-1.0, 1.0) * sampleRadius;
							tp += Fx::randomSampler(-1.0, 1.0) * sampleRadius;
						}

						Vector3f v;
						Fx::MakeCubeVectorFromFaceST(face, sp, tp, &v.x, &v.y, &v.z);
						v = v.unit();

						Color4f sampleColor;

						float theta;
						float gamma;
						hwSunPbsky->computeThetaGamma(v.x, -v.z, v.y, &theta, &gamma);

						if (theta >= 0.0f) {
							hwSunPbsky->computeSunRadiance4(theta, gamma, sampleColor);

							sampleColor.r = Fx::flterrzero(sampleColor.r * sampleScale);
							sampleColor.g = Fx::flterrzero(sampleColor.g * sampleScale);
							sampleColor.b = Fx::flterrzero(sampleColor.b * sampleScale);

							//float r = hwSunPbsky->ComputeSunRadiance2(theta, gamma, 0) * sampleScale;
							//float g = hwSunPbsky->ComputeSunRadiance2(theta, gamma, 1) * sampleScale;
							//float b = hwSunPbsky->ComputeSunRadiance2(theta, gamma, 2) * sampleScale;

							//if (std::finite(r) && std::isfinite(g) && std::isfinite(b))
							//{
							//	sampleColor.r = r;
							//	sampleColor.g = g;
							//	sampleColor.b = b;
							//}
							//else
							//{
							//	sampleColor.r = 0.0f;
							//	sampleColor.g = 0.0f;
							//	sampleColor.b = 0.0f;
							//}
						}
						else {
							sampleColor = groundRadiance * sampleScale;
						}
						// debug to see if cube map is mapped properly
						//sampleColor.r = 0.5f*ptr.X + 0.5f;
						//sampleColor.g = 0.5f*ptr.y + 0.5f;
						//sampleColor.b = 0.5f*ptr.z + 0.5f;
						//sampleColor.r = sp;
						//sampleColor.g = tp;
						//sampleColor.b = 0.0f;
						color += sampleColor;
					}

					if (nSamples > 1)
						color = color / (float)nSamples;

					// flip the t axis
					if (!flipY)
						generatedSunCubeMap.setPixelUnsafe(is, (resolution - 1) - it, face, color);
					else
						generatedSunCubeMap.setPixelUnsafe(is, it, face, color);
				}
			}
		}

		if (normalize) {
			//float average = hwSunPbsky->totalValue / hwSunPbsky->nSamples;
			float invScale = 1.0f / (sampleScale * hwSunPbsky->getMaxValue());
			//generatedSunCubeMap.scaleColors(invScale);
			invScale = 1.0f / hwSunPbsky->getSunDiskRadiance().ToVector3().length();
			generatedSunCubeMap.scaleColors(invScale);
		}
		else {
			// float average = hwSunPbsky->totalValue / hwSunPbsky->nSamples;
			// float invScale = 1.0f / (sampleScale * hwSunPbsky->maxValue);
			// generatedSunCubeMap.scaleColors(invScale);
			// invScale = 1.0f / hwSunPbsky->GetSunDiskRadiance().ToVector3().length();
			// generatedSunCubeMap.scaleColors(2.5f * powf(2.0f, -6.0f));
		}

		minRgbValue = hwSunPbsky->getMinValue();
		maxRgbValue = hwSunPbsky->getMaxValue();
	}

	void PhysicallyBasedSky::computeCylinderMap(int width, int height, bool normalize, float sampleScale) {
		_prepareForCompute();
		generatedSunCylMap.resize(width, height);

		float sampleRadius = nSamples > 1 ? 0.5f : 0.0f;

		// (u, ptr) is the texture_ coordinate in the range [0.0, 1.0]
		// this maps to 0 to 2 * pi and -pi to pi
		// (i, j) are the indices into the texture_ map
		// (theta, gamma) are the coordinates in the sky with (theta, phi) being (theta, pi/2 - gamma)
		//
		int i, j;
		float du = 2.0f / (generatedSunCylMap.width() - 1.0f); // subtract 1.0 from width to ensure image covers -1.0 to 1.0
		float dv = 2.0f / (generatedSunCylMap.height() - 1.0f);

		float u = -1.0f;
		for (i = 0; i < width; i++) {
			float v = -1.0f;
			for (j = 0; j < height; j++) {
				Color4f color;

				for (int s = 0; s < nSamples; s++) {
					Color4f sampleColor;
					float us = (float)(u)+Fx::randomSampler(-du, du) * sampleRadius;
					float vs = (float)(v)+Fx::randomSampler(-dv, dv) * sampleRadius;

					us = std::clamp<float>(us, -1.0f, 1.0f);
					vs = std::clamp<float>(vs, -1.0f, 1.0f);

					float inclination = ((1.0f + vs) / 2.0f * (float)FX_PIOVERTWO);
					float azimuth = (1.0f + us) / 2.0f * (float)FX_TWOPI;

					float theta;
					float gamma;
					hwSunPbsky->computeThetaGamma(inclination, azimuth, &theta, &gamma);

					if (theta >= 0.0f) {
						//sampleColor.r = hwSunPbsky->ComputeSunRadiance2(theta, gamma, 0) * sampleScale;
						//sampleColor.g = hwSunPbsky->ComputeSunRadiance2(theta, gamma, 1) * sampleScale;
						//sampleColor.b = hwSunPbsky->ComputeSunRadiance2(theta, gamma, 2) * sampleScale;

						hwSunPbsky->computeSunRadiance4(theta, gamma, sampleColor);

						sampleColor.r = Fx::flterrzero(sampleColor.r * sampleScale);
						sampleColor.g = Fx::flterrzero(sampleColor.g * sampleScale);
						sampleColor.b = Fx::flterrzero(sampleColor.b * sampleScale);
					}
					color += sampleColor;
				}

				if (nSamples > 1)
					color = color / (float)nSamples;

				generatedSunCylMap.setPixel(i, j, color);
				v += dv;
			}
			u += du;
		}

		if (normalize) {
			float invScale = 1.0f / hwSunPbsky->getMaxValue();
			for (i = 0; i < width; i++) {
				for (j = 0; j < height; j++) {
					Color4f color = generatedSunCylMap.getPixel(i, j);
					color *= invScale;
				}
			}
		}

		minRgbValue = hwSunPbsky->getMinValue();
		maxRgbValue = hwSunPbsky->getMaxValue();
	}


	//void PhysicallyBasedSky::ComputeSphereMap(int width, int height, bool normalize, float sampleScale)
	//{
	//}


	void PhysicallyBasedSky::computeSunGroundRadiances() {
		_prepareForCompute(false);
		sunDiskRadiance = hwSunPbsky->getSunDiskRadiance();
		groundRadiance = hwSunPbsky->getGroundRadiance();
	}


	void PhysicallyBasedSky::_prepareForCompute(bool resetStats) {
		minRgbValue = FLT_MAX;
		maxRgbValue = -FLT_MAX;
		hwSunPbsky->init(turbidity, groundAlbedo, static_cast<float>(sunPosition_.a), static_cast<float>(sunPosition_.A));
		if (resetStats)
			hwSunPbsky->resetStatisticSamples();
		hwMoonPbsky->init(turbidity, groundAlbedo, (float)moonPosition_.a, (float)moonPosition_.A);
	}


	Color3f PhysicallyBasedSky::computeModisAlbedo(bool recalc) const {
		if (recalc || recalcModis_) {
			computeModisAlbedo((float)astroCalc.getLatitude(), (float)astroCalc.getLongitude(), (float)astroCalc.getMonthOfYear());
		}
		return modisColor_;
	}


	Color3f PhysicallyBasedSky::computeModisAlbedo(float latitude, float longitude, float month) const {
		constexpr bool hires = true;
		using string_vector = std::vector<std::string>;
		static string_vector modis_data_294x196{
			"ssphh-data/resources/textures/world.200401x294x196.jpg",
			"ssphh-data/resources/textures/world.200402x294x196.jpg",
			"ssphh-data/resources/textures/world.200403x294x196.jpg",
			"ssphh-data/resources/textures/world.200404x294x196.jpg",
			"ssphh-data/resources/textures/world.200405x294x196.jpg",
			"ssphh-data/resources/textures/world.200406x294x196.jpg",
			"ssphh-data/resources/textures/world.200407x294x196.jpg",
			"ssphh-data/resources/textures/world.200408x294x196.jpg",
			"ssphh-data/resources/textures/world.200409x294x196.jpg",
			"ssphh-data/resources/textures/world.200410x294x196.jpg",
			"ssphh-data/resources/textures/world.200411x294x196.jpg",
			"ssphh-data/resources/textures/world.200412x294x196.jpg",
		};
		static string_vector modis_data_21600x10800{
			"ssphh-data/resources/textures/world.200401.3x21600x10800.jpg",
			"ssphh-data/resources/textures/world.200402.3x21600x10800.jpg",
			"ssphh-data/resources/textures/world.200403.3x21600x10800.jpg",
			"ssphh-data/resources/textures/world.200404.3x21600x10800.jpg",
			"ssphh-data/resources/textures/world.200406.3x21600x10800.jpg",
			"ssphh-data/resources/textures/world.200407.3x21600x10800.jpg",
			"ssphh-data/resources/textures/world.200408.3x21600x10800.jpg",
			"ssphh-data/resources/textures/world.200409.3x21600x10800.jpg",
			"ssphh-data/resources/textures/world.200410.3x21600x10800.jpg",
			"ssphh-data/resources/textures/world.200411.3x21600x10800.jpg",
			"ssphh-data/resources/textures/world.200412.3x21600x10800.jpg",
		};
		static constexpr int width = hires ? 21600 : 294;
		static constexpr int height = hires ? 10800 : 196;
		static constexpr int width_over_2 = width >> 1;
		static constexpr int height_over_2 = height >> 1;
		static std::vector<SDL_Surface*> surfaces(12, nullptr);
		static std::vector<Image3f> modis(12);
		float fract = month - std::floor(month);
		int month1 = (int)(month - fract - 1);
		int month2 = month1 + 1;
		month1 = month1 % 12;
		month2 = month2 % 12;
		latitude = std::clamp(latitude, -90.0f, 90.0f) / 180.0f + 0.5f;
		longitude = Fx::wrap(longitude, -180.0f, 180.0f) / 360.0f + 0.5f;
		int s = std::clamp((int)((height - .001f) * latitude), 0, height - 1);
		int t = std::clamp((int)((width - 0.001f) * longitude), 0, width - 1);
		Image3f& s1 = modis[month1];
		Image3f& s2 = modis[month2];
		if (!s1) {
			Hf::StopWatch stopwatch;
			const char* file = hires ? modis_data_21600x10800[month1].c_str() : modis_data_294x196[month1].c_str();
			if (hires)
				LoadImage3f(file, s1);
			else
				LoadImage3f(file, s1);
			HFLOGDEBUG("Loaded %s in %3.2fms", file, stopwatch.Stop_msf());
		}
		if (!s2) {
			Hf::StopWatch stopwatch;
			const char* file = hires ? modis_data_21600x10800[month2].c_str() : modis_data_294x196[month2].c_str();
			if (hires)
				LoadImage3f(file, s2);
			else
				LoadImage3f(file, s2);
			HFLOGDEBUG("Loaded %s in %3.2fms", file, stopwatch.Stop_msf());
		}
		if (!s1 || !s2) return { 0.0f, 0.0f, 0.0f };
		Color3f p1 = s1.getPixel(s, t);
		Color3f p2 = s2.getPixel(s, t);
		modisColor_ = lerp(fract, p1, p2);
		recalcModis_ = false;
		return modisColor_;
	}

} // namespace Sf

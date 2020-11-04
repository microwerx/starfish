#ifndef FLUXIONS_PBSKY_HPP
#define FLUXIONS_PBSKY_HPP

#include <starfish-astronomy.hpp>
#include <starfish-hosek-wilkie.hpp>
#include <fluxions_gte.hpp>
#include <fluxions_gte_image.hpp>
#include <starfish-ipbsky.hpp>

namespace Sf {
	namespace Fx = ::Fluxions;

	/////////////////////////////////////////////////////////////////////
	// P H Y S I C A L L Y   B A S E D   S K Y //////////////////////////
	/////////////////////////////////////////////////////////////////////

	//class HosekWilkiePBSky;


	class PhysicallyBasedSky {
	public:
		PhysicallyBasedSky();
		~PhysicallyBasedSky();

		void setLocation(float latitude, float longitude);
		float getLatitude() const { return (float)astroCalc.getLatitude(); }
		float getLongitude() const { return (float)astroCalc.getLongitude(); }
		const PA::CivilDateTime getCivilDateTime() const { return astroCalc.GetDateTime(); }
		time_t getTime() const;
		void setCivilDateTime(const PA::CivilDateTime& dtg);
		void setTime(time_t t, float fractSeconds = 0.0f);
		void setLocalDate(int day, int month, int year, bool isdst, int timeOffset);
		void setLocalTime(int hh, int mm, int ss, float ss_frac);
		void setTurbidity(float T);
		float getTurbidity() const;
		void setSunPosition(double azimuth, double altitude);
		void setSunPosition(double sunLong);
		void setMoonPosition(double RA, double dec);
		void setMoonPosition(EquatorialCoord moonRADec);
		float getAverageRadiance() const;
		float getSunAzimuth() const { return (float)sunPosition_.A; }
		float getSunAltitude() const { return (float)sunPosition_.a; }
		void setGroundAlbedo(float r, float g, float b) { groundAlbedo.reset(r, g, b); };
		const Fx::Color4f& getGroundAlbedo() const { return groundAlbedo; }

		// Compute MODIS Albedo at the specified latitude and longitude, measured in degrees
		Fx::Color3f computeModisAlbedo(bool recalc) const;
		Fx::Color3f computeModisAlbedo(float latitude, float longitude, float month) const;

		float getDayOfYear() const { return (float)astroCalc.getDayOfYear(); }
		float getMonthOfYear() const { return (float)astroCalc.getMonthOfYear(); }

		void setNumSamples(int samples) { nSamples = std::clamp<int>(samples, 1, 16); }
		int getNumSamples() const { return nSamples; }

		void computeAstroFromLocale();

		void computeCubeMap(int resolution, bool normalize = false, float sampleScale = 8.0f, bool flipY = false);
		// not implemented
		void computeCylinderMap(int width, int height, bool normalize = false, float sampleScale = 8.0f);
		// not implemented
		//void computeSphereMap(int width, int height, bool normalize = false, float sampleScale = 8.0f);

		void computeSunGroundRadiances();
		Fx::Color4f getSunDiskRadiance() const { return sunDiskRadiance; }
		Fx::Color4f getGroundRadiance() const { return groundRadiance; }

		float getMinRgbValue() const { return minRgbValue; }
		float getMaxRgbValue() const { return maxRgbValue; }

		Fx::Image4f generatedSunCubeMap;
		// unused
		Fx::Image4f generatedSunSphMap;
		// unused
		Fx::Image4f generatedSunCylMap;

		Fx::Image4f generatedMoonCubeMap;
		Fx::Image4f generagedMoonCylMap;

		int getDay() const { return astroCalc.GetDateTime().day; }
		int getMonth() const { return astroCalc.GetDateTime().month; }
		int getYear() const { return astroCalc.GetDateTime().year; }
		int getHour() const { return astroCalc.GetDateTime().hh; }
		int getMin() const { return astroCalc.GetDateTime().mm; }
		int getSec() const { return astroCalc.GetDateTime().ss; }
		double getSecFract() const { return astroCalc.GetDateTime().ss_frac; }
		double getLST() const { return astroCalc.getLST(); }

		Fx::Vector3f sunDirTo() const { return  sunVector_; }
		Fx::Vector3f moonDirTo() const { return moonVector_; }
		EquatorialCoord getSunRADec() const { return sunRADec_; }
		EquatorialCoord getMoonRADec() const { return moonRADec_; }
		HorizonCoord getSunPosition() const { return sunPosition_; }
		HorizonCoord getMoonPosition() const { return moonPosition_; }
	private:
		void _computeSunFromLocale();
		void _computeMoonFromLocale();

		mutable bool recalcModis_{ true };
		mutable Fx::Color3f modisColor_;

		AstroCalc astroCalc;
		HorizonCoord sunPosition_;
		HorizonCoord moonPosition_;
		EquatorialCoord sunRADec_;
		EquatorialCoord moonRADec_;
		Fx::Vector3f sunVector_;
		Fx::Vector3f moonVector_;
		float turbidity = 1.0f;
		Fx::Color4f groundAlbedo;
		Fx::Color4f groundRadiance;
		Fx::Color4f sunDiskRadiance;
		int nSamples;
		float minRgbValue;
		float maxRgbValue;
		//std::shared_ptr<HosekWilkiePBSky> hwSunPbsky;
		//std::shared_ptr<HosekWilkiePBSky> hwMoonPbsky;
		std::shared_ptr<IPBSky> hwSunPbsky;
		std::shared_ptr<IPBSky> hwMoonPbsky;

		void _prepareForCompute(bool resetStats = true);
	};
} // namespace Fluxions

#endif

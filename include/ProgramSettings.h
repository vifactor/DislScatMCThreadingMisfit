/*
 * ProgramSettings.h
 *
 *  Created on: 11 june 2013
 *      Author: kopp
 */

#ifndef PROGRAMSETTINGS_H_
#define PROGRAMSETTINGS_H_

#include "StringTools.h"
#include "MillerIndexHex.h"
#include <libconfig.h++>

struct Range
{
	double m_min;
	double m_max;
	size_t m_sampling;

	double getStep() const
	{
		return (m_sampling > 1) ? (m_max - m_min) / (m_sampling - 1) : 0.0;
	}
	void toVector(std::vector<double>& vec) const
	{
		double step = getStep();
		for (size_t i = 0; i < m_sampling; ++i)
		{
			vec.push_back(m_min + step * i);
		}
	}
	bool good() const
	{
		if ((m_min <= m_max) && (m_sampling > 0))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

class ProgramSettings
{
public:
	class Exception: public std::exception
	{
	public:
		Exception(std::string m)
		{
			msg = "ProgramSettings::" + m;
		}
		~Exception() throw ()
		{
		}
		const char* what() const throw ()
		{
			return msg.c_str();
		}
	private:
		std::string msg;
	};
	struct SampleSettings
	{
		double nu;
		double thickness;
		double width;
		bool isHalfSpace;
		/*hexagonal lattice parameters*/
		double a0, c0;

		struct MisfitDislocationInterface
		{
			/*burgers components*/
			double b_x, b_z;
			double rho;
		};
        struct MisfitDislocationFamily
		{
		    /*dislocation density*/
			double rho;
			/*burgers vectors*/
			std::vector<MillerDirectHexIndices> bs;
			/*dislocation lines*/
			std::vector<MillerDirectHexIndices> ls;
			/*order parameter*/
			int gamma;
		};
		struct ThreadingDislocationFamily
		{
			std::vector<MillerDirectHexIndices> bs;
			double rho;
			double rc;
		};

		std::vector<MisfitDislocationInterface> misfit_interface;
		std::vector<MisfitDislocationFamily> misfit_family;
		std::vector<ThreadingDislocationFamily> threading_family;
	};
	struct CalculatorSettings
	{
		/*hexagonal index*/
		MillerReciprocalHexIndices Q;
		/*X-ray wavelength*/
		double lambda;

		double qresolX, qresolZ;
	};
	struct EngineSettings
	{
		std::string outfile;
		Range qxRange, qzRange;
		unsigned long long nbMCCalls;
        double precision;

		enum Geometry {geomSKEW, geomCOPLANAR, geomUNKNOWN} m_geometry;
		enum Diffractometry {diffDOUBLE, diffTRIPLE, diffUNKNOWN} m_diffractometry;
	};
	const SampleSettings& getSampleSettings() const
	{
		return m_sampleSettings;
	}
	const CalculatorSettings& getCalculatorSettings() const
	{
		return m_calculatorSettings;
	}
	const EngineSettings& getEngineSettings() const
	{
		return m_engineSettings;
	}
	const std::string& getConfigfile() const
	{
		return m_cfgfile;
	}
	ProgramSettings();
	virtual ~ProgramSettings();

	void read(const std::string& cfg);
	void print() const;
protected:
	void readCalculatorSettings(const libconfig::Setting& root);
	void readSampleSettings(const libconfig::Setting& root);
	void readEngineSettings(const libconfig::Setting& root);

	void readCoplanarSettings(const libconfig::Setting& stg);
	void readSkewSettings(const libconfig::Setting& stg);
	void readSkewDoubleSettings(const libconfig::Setting& stg);
	void readSkewTripleSettings(const libconfig::Setting& stg);

	void readMisfitDislocationInterfaces(const libconfig::Setting& stg);
	void readThreadingDislocationFamilies(const libconfig::Setting& stg);
	void readMisfitDislocationFamilies(const libconfig::Setting& stg);

	EngineSettings::Diffractometry defineDiffractometry(const libconfig::Setting& stg);
	EngineSettings::Geometry defineGeometry(const libconfig::Setting& stg);

	void printSampleSettings() const;
	void printCalculatorSettings() const;
	void printEngineSettings() const;

	void printCoplanarSettings() const;
	void printSkewSettings() const;
	void printSkewDoubleSettings() const;
	void printSkewTripleSettings() const;

	void printMisfitDislocationInterfaces() const;
	void printThreadingDislocationFamilies() const;
	void printMisfitDislocationFamilies() const;

	SampleSettings m_sampleSettings;
	CalculatorSettings m_calculatorSettings;
	EngineSettings m_engineSettings;
	std::string m_cfgfile;
};

#endif /* PROGRAMSETTINGS_H_ */

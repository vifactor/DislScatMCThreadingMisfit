/*
 * ProgramSettings.cpp
 *
 *  Created on: 11 june 2013
 *      Author: kopp
 */

#include "ProgramSettings.h"

Range readRange(const libconfig::Setting& stg)
{
	Range range;

	range.m_min = stg[0][0];
	range.m_max = stg[0][1];
	range.m_sampling = stg[1];

	return range;
}

std::ostream& operator<<(std::ostream& out, const Range& range)
{
	out << "[" << range.m_min << ", " << range.m_max << "]:" << range.m_sampling;
	return out;
}

MillerDirectHexIndices readMillerDirectHexIndices(const libconfig::Setting& stg)
{
    MillerDirectHexIndices index;
	if(stg.isArray() && stg.getLength() == MillerHexIndicesDimension)
	{
			index.X = stg[0];
			index.Y = stg[1];
			index.T = stg[2];
			index.Z = stg[3];
	}
	else
	{
		throw ProgramSettings::Exception("Check setting: " + toString(stg.getPath()));
	}
	return index;
}

MillerReciprocalHexIndices readMillerReciprocalHexIndices(const libconfig::Setting& stg)
{
    MillerReciprocalHexIndices index;
	if(stg.isArray() && stg.getLength() == MillerHexIndicesDimension)
	{
			index.H = stg[0];
			index.K = stg[1];
			index.I = stg[2];
			index.L = stg[3];
	}
	else
	{
		throw ProgramSettings::Exception("Check setting: " + toString(stg.getPath()));
	}
	return index;
}

std::ostream& operator<<(std::ostream& out, const MillerDirectHexIndices& index)
{
	out << "[" << index.X << ", " << index.Y << ", " << index.T << ", " << index.Z << "]";
	return out;
}

std::ostream& operator<<(std::ostream& out, const MillerReciprocalHexIndices& index)
{
	out << "[" << index.H << ", " << index.K << ", " << index.I << ", " << index.L << "]";
	return out;
}

ProgramSettings::ProgramSettings()
{
}

ProgramSettings::~ProgramSettings()
{
}

void ProgramSettings::read(const std::string& cfgfile)
{
	libconfig::Config cfg;

	m_cfgfile = cfgfile;
	// Read the file. If there is an error, report it
	try
	{
		cfg.readFile(m_cfgfile.c_str());
		cfg.setAutoConvert(true);
		const libconfig::Setting& root = cfg.getRoot();

		readSampleSettings(root);
		readCalculatorSettings(root);
		readEngineSettings(root);
	} catch (const libconfig::FileIOException &fioex)
	{
		throw Exception(toString(fioex.what()) + " in\t" + cfgfile);
	} catch (const libconfig::ParseException &pex)
	{
		throw Exception(
				toString(pex.what()) + " in\t" + cfgfile + ":"
						+ toString(pex.getLine()) + " - "
						+ toString(pex.getError()));
	} catch (const libconfig::SettingNotFoundException &nfex)
	{
		throw Exception(
				toString(nfex.what()) + "\t" + toString(nfex.getPath())
						+ " in\t" + cfgfile);
	} catch (libconfig::SettingTypeException& tex)
	{
		throw Exception(
				toString(tex.what()) + "\t" + toString(tex.getPath()) + " in\t"
						+ cfgfile);
	}
}

void ProgramSettings::readCalculatorSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &calculator = root["Calculator"];

	/*reflection*/
	if(calculator["Q"].isArray() && calculator["Q"].getLength() == MillerHexIndicesDimension)
	{
			m_calculatorSettings.Q.H = calculator["Q"][0];
			m_calculatorSettings.Q.K = calculator["Q"][1];
			m_calculatorSettings.Q.I = calculator["Q"][2];
			m_calculatorSettings.Q.L = calculator["Q"][3];
	}
	else
	{
		throw ProgramSettings::Exception("Check reflection: " + toString(calculator["Q"].getPath()));
	}
	/*check the property of hexagonal Miller indices*/
	if(!m_calculatorSettings.Q.isCorrect())
	{
		throw ProgramSettings::Exception( "Miller index property is brocken: " + toString(calculator["Q"].getPath()));
	}

	/*X-ray wavelength*/
	m_calculatorSettings.lambda = calculator["lambda"];

	/*m_calculatorSettings.scale = calculator["scale"];
	m_calculatorSettings.background = calculator["background"];*/

	m_calculatorSettings.qresolX = calculator["resolution"]["x"];
	m_calculatorSettings.qresolZ = calculator["resolution"]["z"];
}

void ProgramSettings::readSampleSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &sample = root["Sample"];

	/*lattice parameters*/
	m_sampleSettings.a0 = sample["a0"];
	m_sampleSettings.c0 = sample["c0"];

	/*Poisson ratio*/
	m_sampleSettings.nu = sample["nu"];

	/*Sample sizes*/
	m_sampleSettings.thickness = sample["thickness"];
	m_sampleSettings.width = sample["width"];

	/*should threading dislocations be considered with surface relaxation term*/
	m_sampleSettings.isHalfSpace = sample["isHalfSpace"];

	/*dislocation settings*/
	const libconfig::Setting &dislocations = sample["dislocations"];

	if(dislocations.exists("misfit_interfaces"))
	{
		const libconfig::Setting &misfit_interfaces = dislocations["misfit_interfaces"];
		readMisfitDislocationInterfaces(misfit_interfaces);
	}

    if(dislocations.exists("misfit_families"))
	{
		const libconfig::Setting &misfit_families = dislocations["misfit_families"];
		readMisfitDislocationFamilies(misfit_families);
	}

	if(dislocations.exists("threading_families"))
	{
		const libconfig::Setting &threading_families = dislocations["threading_families"];
		readThreadingDislocationFamilies(threading_families);
	}
}

void ProgramSettings::readEngineSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &engine = root["Engine"];

	m_engineSettings.m_geometry = defineGeometry(engine["geometry"]);
	m_engineSettings.nbMCCalls = engine["nbMCCalls"];
	m_engineSettings.precision = engine["precision"];

	if(m_engineSettings.m_geometry == EngineSettings::geomCOPLANAR)
	{
		const libconfig::Setting &settings = engine["coplanar_settings"];
		readCoplanarSettings(settings);
	}
	else if(m_engineSettings.m_geometry == EngineSettings::geomSKEW)
	{
		const libconfig::Setting &settings = engine["skew_settings"];
		readSkewSettings(settings);
	}
	else
	{
		throw Exception("Unknown geometry setup:\t" + toString(engine["geometry"].c_str()));
	}

	m_engineSettings.outfile = engine["outfile"].c_str();
}

void ProgramSettings::readCoplanarSettings(const libconfig::Setting& stg)
{
	m_engineSettings.qxRange = readRange(stg["qxrange"]);
	m_engineSettings.qzRange = readRange(stg["qzrange"]);
}

void ProgramSettings::readSkewSettings(const libconfig::Setting& stg)
{
	m_engineSettings.m_diffractometry = defineDiffractometry(stg["diffractometry"]);

	if(m_engineSettings.m_diffractometry == EngineSettings::diffDOUBLE)
	{
		readSkewDoubleSettings(stg["double_settings"]);
	}else if(m_engineSettings.m_diffractometry == EngineSettings::diffTRIPLE)
	{
		readSkewTripleSettings(stg["triple_settings"]);
	}
	else
	{
		throw Exception("Unknown diffractometry setup:\t" + toString(stg["diffractometry"].c_str()));
	}
}

ProgramSettings::EngineSettings::Diffractometry ProgramSettings::defineDiffractometry(const libconfig::Setting& stg)
{
	std::string diffractometry;

	diffractometry = stg.c_str();
	if (diffractometry.compare("DOUBLE") == 0)
	{
		return EngineSettings::diffDOUBLE;
	}
	else if (diffractometry.compare("TRIPLE") == 0)
	{
		return EngineSettings::diffTRIPLE;
	}
	else
	{
		return EngineSettings::diffUNKNOWN;
	}
}

ProgramSettings::EngineSettings::Geometry ProgramSettings::defineGeometry(const libconfig::Setting& stg)
{
	std::string diffractometry;

	diffractometry = stg.c_str();
	if (diffractometry.compare("COPLANAR") == 0)
	{
		return EngineSettings::geomCOPLANAR;
	}
	else if (diffractometry.compare("SKEW") == 0)
	{
		return EngineSettings::geomSKEW;
	}
	else
	{
		return EngineSettings::geomUNKNOWN;
	}
}

void ProgramSettings::print() const
{
	printEngineSettings();
	printSampleSettings();
	printCalculatorSettings();
}

void ProgramSettings::printCalculatorSettings() const
{
	std::cout << "---Calculator settings---" << std::endl;
	std::cout << "Reflection:\t" << m_calculatorSettings.Q << std::endl;
	std::cout << "X-ray wavelength:\t" << m_calculatorSettings.lambda << std::endl;
	std::cout << "Resolutions (dqx, dqz):\t" << m_calculatorSettings.qresolX
			<< "\t" << m_calculatorSettings.qresolZ << std::endl;
	//std::cout << "Intensity scale coefficient:\t" << m_calculatorSettings.scale << std::endl;
	//std::cout << "Intensity background:\t" << m_calculatorSettings.background << std::endl;
}

void ProgramSettings::printSampleSettings() const
{
	std::cout << "---Sample settings---" << std::endl;
	std::cout << "Lattice parameters: (a0, c0)\t" << m_sampleSettings.a0 << ", "
			<< m_sampleSettings.c0 << std::endl;
	std::cout << "Sample sizes (thickness width):\t" << m_sampleSettings.thickness << "\t"
			<< m_sampleSettings.width << std::endl;
	std::cout << "Poisson ratio:\t" << m_sampleSettings.nu << std::endl;

	std::cout << "Misfit dislocation interfaces:" << std::endl;
	printMisfitDislocationInterfaces();

	std::cout << "Threading dislocation families:" << std::endl;
	printThreadingDislocationFamilies();

    std::cout << "Misfit dislocation families:" << std::endl;
	printMisfitDislocationFamilies();

}

void ProgramSettings::printEngineSettings() const
{
	std::cout << "---Engine settings---" << std::endl;

	switch(m_engineSettings.m_geometry)
	{
	case EngineSettings::geomCOPLANAR:
		std::cout << "Intensity is calculated in co-planar geometry" << std::endl;
		printCoplanarSettings();
		break;
	case EngineSettings::geomSKEW:
		std::cout << "Intensity is calculated in skew geometry" << std::endl;
		printSkewSettings();
		break;
	default:
		//never happens
		break;
	}

    std::cout << "Stop criteria :\t" << m_engineSettings.precision << std::endl;
    std::cout << "Initial number of MC steps before output:\t" << m_engineSettings.nbMCCalls << std::endl;
	std::cout << "Output basename:\t" << m_engineSettings.outfile << std::endl;
}

void ProgramSettings::printCoplanarSettings() const
{
	std::cout << "Qx range:\t" << m_engineSettings.qxRange << std::endl;
	std::cout << "Qz range:\t" << m_engineSettings.qzRange << std::endl;
	std::cout << "Initial nb Monte-Carlo steps:\t" << m_engineSettings.nbMCCalls << std::endl;
}

void ProgramSettings::printSkewSettings() const
{
	switch(m_engineSettings.m_diffractometry)
	{
	case EngineSettings::diffDOUBLE:
		std::cout << "Double crystal diffractometry is used" << std::endl;
		std::cout << "Q range:\t" << m_engineSettings.qzRange << std::endl;
		break;
	case EngineSettings::diffTRIPLE:
		std::cout << "Triple crystal diffractometry is used" << std::endl;
		std::cout << "Qx range:\t" << m_engineSettings.qxRange << std::endl;
		std::cout << "Qz range:\t" << m_engineSettings.qzRange << std::endl;
		break;
	default:
		//never happens
		break;
	}
}

void ProgramSettings::readSkewDoubleSettings(const libconfig::Setting& stg)
{
	m_engineSettings.qzRange = readRange(stg["qrange"]);
}

void ProgramSettings::readSkewTripleSettings(const libconfig::Setting& stg)
{
	m_engineSettings.qxRange = readRange(stg["qxrange"]);
	m_engineSettings.qzRange = readRange(stg["qzrange"]);
}

void ProgramSettings::readMisfitDislocationInterfaces(const libconfig::Setting& stg)
{
	size_t nbDislocations;

	if(stg.isList())
	{
		nbDislocations = stg.getLength();
		m_sampleSettings.misfit_interface.resize(nbDislocations);
		for(size_t i = 0; i < nbDislocations; ++i)
		{
			m_sampleSettings.misfit_interface[i].rho = stg[i]["rho"];
			m_sampleSettings.misfit_interface[i].b_x = stg[i]["b_x"];
			m_sampleSettings.misfit_interface[i].b_z = stg[i]["b_z"];
		}
	}
}

void ProgramSettings::readMisfitDislocationFamilies(const libconfig::Setting& stg)
{
	size_t nb_families;
	size_t nb_vectors;
	MillerDirectHexIndices b, l;

	if(stg.isList())
	{
		nb_families = stg.getLength();
		m_sampleSettings.misfit_family.resize(nb_families);
		for(size_t ifam = 0; ifam < nb_families; ++ifam)
		{
            m_sampleSettings.misfit_family[ifam].rho = stg[ifam]["rho"];
            m_sampleSettings.misfit_family[ifam].gamma = stg[ifam]["gamma"];
            const libconfig::Setting& vecs = stg[ifam]["vectors"];
            nb_vectors = vecs.getLength();
            for(size_t ivec = 0; ivec < nb_vectors; ++ivec)
            {
                b = readMillerDirectHexIndices(vecs[ivec]["b"]);
                l = readMillerDirectHexIndices(vecs[ivec]["l"]);

                /*check the property of hexagonal Miller indices*/
                if(!b.isCorrect())
                {
                    throw ProgramSettings::Exception("Miller hexagonal index: " + toString(vecs[ivec]["b"].getPath()));
                }
                if(!l.isCorrect())
                {
                    throw ProgramSettings::Exception("Miller hexagonal index: " + toString(vecs[ivec]["l"].getPath()));
                }

                m_sampleSettings.misfit_family[ifam].bs.push_back(b);
                m_sampleSettings.misfit_family[ifam].ls.push_back(l);
            }
		}
	}
}

void ProgramSettings::readThreadingDislocationFamilies(const libconfig::Setting& stg)
{
	size_t nb_families;
	size_t nb_vectors;
	MillerDirectHexIndices b;

	if(stg.isList())
	{
		nb_families = stg.getLength();
		m_sampleSettings.threading_family.resize(nb_families);
		for(size_t ifam = 0; ifam < nb_families; ++ifam)
		{
            m_sampleSettings.threading_family[ifam].rho = stg[ifam]["rho"];
            m_sampleSettings.threading_family[ifam].rc = stg[ifam]["rc"];
            const libconfig::Setting& vecs = stg[ifam]["vectors"];
            nb_vectors = vecs.getLength();
            for(size_t ivec = 0; ivec < nb_vectors; ++ivec)
            {
                b = readMillerDirectHexIndices(vecs[ivec]);

                /*check the property of hexagonal Miller indices*/
                if(!b.isCorrect())
                {
                    throw ProgramSettings::Exception("Miller hexagonal index: " + toString(vecs[ivec].getPath()));
                }

                m_sampleSettings.threading_family[ifam].bs.push_back(b);
            }
		}
	}
}

void ProgramSettings::printMisfitDislocationInterfaces() const
{
	for (size_t i = 0; i < m_sampleSettings.misfit_interface.size(); ++i)
	{
		std::cout << "\tSet:\t" << i << std::endl;
		std::cout << "\t\tBurgers vector (bx):\t"
				<< m_sampleSettings.misfit_interface[i].b_x << std::endl;
		std::cout << "\t\tBurgers vector (bz):\t"
				<< m_sampleSettings.misfit_interface[i].b_z << std::endl;
		std::cout << "\t\tDensity :\t" << m_sampleSettings.misfit_interface[i].rho
				<< std::endl;
	}
}

void ProgramSettings::printThreadingDislocationFamilies() const
{
	for (size_t ifam = 0; ifam < m_sampleSettings.threading_family.size(); ++ifam)
	{
		std::cout << "\tFamily:\t" << ifam << std::endl;
		std::cout << "\t\tDensity:\t" << m_sampleSettings.threading_family[ifam].rho
				<< std::endl;
		std::cout << "\t\tCorrelation radius:\t" << m_sampleSettings.threading_family[ifam].rc
				<< std::endl;

        std::cout << "\t\tBurgers vectors:\t" << std::endl;
        for (size_t ivec = 0; ivec < m_sampleSettings.threading_family[ifam].bs.size(); ++ivec)
        {
             std::cout << "\t\t\t" << m_sampleSettings.threading_family[ifam].bs[ivec] << std::endl;
        }
	}
}

void ProgramSettings::printMisfitDislocationFamilies() const
{
	for (size_t ifam = 0; ifam < m_sampleSettings.misfit_family.size(); ++ifam)
	{
		std::cout << "\tFamily:\t" << ifam << std::endl;
		std::cout << "\t\tDensity:\t" << m_sampleSettings.misfit_family[ifam].rho
				<< std::endl;
		std::cout << "\t\tCorrelation parameter:\t" << m_sampleSettings.misfit_family[ifam].gamma
				<< std::endl;

        std::cout << "\t\tBurgers vectors:Dislocation lines:\t" << std::endl;
        for (size_t ivec = 0; ivec < m_sampleSettings.misfit_family[ifam].bs.size(); ++ivec)
        {
             std::cout << "\t\t\tb" << m_sampleSettings.misfit_family[ifam].bs[ivec] << ":l" <<
                    m_sampleSettings.misfit_family[ifam].ls[ivec] << std::endl;
        }
	}
}

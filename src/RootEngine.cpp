#include "RootEngine.h"

RootEngine::RootEngine()
{
    m_scale = 1.0;
}

RootEngine::~RootEngine()
{
}

void RootEngine::setup(const std::string& cfgfile, int seed)
{
    Engine::setup(cfgfile, seed);
    copySettings(cfgfile, m_programSettings.getEngineSettings().outfile);
    m_programSettings.print();
    std::cout << "Total nb dislocations:\t" << m_sample->getNbDislocations() << std::endl;
    std::cout << "Nb misfit dislocations :\t" << m_sample->getNbMfDislocations() << std::endl;
    std::cout << "Nb threading dislocations :\t" << m_sample->getNbThDislocations() << std::endl;

    switch(m_programSettings.getEngineSettings().m_geometry)
    {
    case ProgramSettings::EngineSettings::geomCOPLANAR:
        setupCoplanarCoefficients();
        break;
    case ProgramSettings::EngineSettings::geomSKEW:
        if(m_programSettings.getEngineSettings().m_diffractometry == ProgramSettings::EngineSettings::diffDOUBLE)
        {
            setupSkewDoubleCoefficients();

        } else if(m_programSettings.getEngineSettings().m_diffractometry == ProgramSettings::EngineSettings::diffTRIPLE)
        {
            setupSkewTripleCoefficients();
        }
        break;
    default:
        break;
    }
}

void RootEngine::setupCoplanarCoefficients()
{

    m_scale = 2 * M_PI * m_sample->m_thickness /
            (m_programSettings.getCalculatorSettings().qresolX *
                m_programSettings.getCalculatorSettings().qresolZ);
}

void RootEngine::setupSkewTripleCoefficients()
{

    m_scale = 2 * M_PI * m_sample->m_thickness /
            (m_programSettings.getCalculatorSettings().qresolX *
                m_programSettings.getCalculatorSettings().qresolZ);
}

void RootEngine::setupSkewDoubleCoefficients()
{
    double sin_theta, Qx, Qz;
    const MillerReciprocalHexIndices& Q = m_programSettings.getCalculatorSettings().Q;

    Qx = 2 * M_PI * sqrt(2.0 / 3 * (Q.H * Q.H + Q.K * Q.K + Q.I * Q.I))
	  / m_programSettings.getSampleSettings().a0;
	Qz = 2 * M_PI * Q.L / m_programSettings.getSampleSettings().c0;

    m_Q_norm = sqrt(Qx * Qx + Qz * Qz);

    sin_theta = m_Q_norm * m_programSettings.getCalculatorSettings().lambda / (4. * M_PI);

    m_cos_theta = sqrt(1.0 - gsl_pow_2(sin_theta));

    m_scale = sqrt(2 * M_PI) * m_sample->m_thickness /
                    m_programSettings.getCalculatorSettings().qresolZ;

}

void RootEngine::save(const MCCalculator::Data * data)
{
    switch(m_programSettings.getEngineSettings().m_geometry)
    {
    case ProgramSettings::EngineSettings::geomCOPLANAR:
        saveMCCalculatorCoplanarData(data);
        break;
    case ProgramSettings::EngineSettings::geomSKEW:
        if(m_programSettings.getEngineSettings().m_diffractometry == ProgramSettings::EngineSettings::diffDOUBLE)
        {
            saveMCCalculatorSkewDoubleData(data);
        } else if(m_programSettings.getEngineSettings().m_diffractometry == ProgramSettings::EngineSettings::diffTRIPLE)
        {
            saveMCCalculatorSkewTripleData(data);
        }
        break;
    default:
        break;
    }
}

void RootEngine::saveMCCalculatorCoplanarData(const MCCalculator::Data * data)
{
    std::ofstream fout;

    fout.open(m_programSettings.getEngineSettings().outfile.c_str());
    fout << "qx\tqz\tre_intensity\tim_intensity" << std::endl;
    fout << "#nbSteps = " << data->m_nbSteps << "\t";
    fout << "ratio = " << data->m_imIMax / data->m_reIMax << std::endl;
    for(size_t ipt = 0; ipt < data->size(); ++ipt)
    {
        fout << data->qx(ipt) << "\t" <<
            data->qz(ipt) << "\t" <<
            m_scale * data->reI(ipt) << "\t" <<
            m_scale * data->imI(ipt) << std::endl;
    }
    fout.close();

    std::cout << "Saved:\t" << data->m_nbSteps << std::endl;
}

void RootEngine::saveMCCalculatorSkewTripleData(const MCCalculator::Data * data)
{
    std::ofstream fout;

    fout.open(m_programSettings.getEngineSettings().outfile.c_str());
    fout << "qx\tqz\tre_intensity\tim_intensity" << std::endl;
    fout << "#nbSteps = " << data->m_nbSteps << "\t";
    fout << "ratio = " << data->m_imIMax / data->m_reIMax << std::endl;
    for(size_t ipt = 0; ipt < data->size(); ++ipt)
    {
        fout << data->qx(ipt) << "\t" <<
            data->qz(ipt) << "\t" <<
            m_scale * data->reI(ipt) << "\t" <<
            m_scale * data->imI(ipt) << std::endl;
    }
    fout.close();

    std::cout << "Saved:\t" << data->m_nbSteps << std::endl;
}

void RootEngine::saveMCCalculatorSkewDoubleData(const MCCalculator::Data * data)
{
    std::ofstream fout;
    double omega;

    fout.open(m_programSettings.getEngineSettings().outfile.c_str());
    fout << "omega\tqz\tre_intensity\tim_intensity" << std::endl;
    fout << "#nbSteps = " << data->m_nbSteps << "\t";
    fout << "ratio = " << data->m_imIMax / data->m_reIMax << std::endl;
    for(size_t ipt = 0; ipt < data->size(); ++ipt)
    {
        omega = (180. / M_PI) * data->qz(ipt)/ ( m_Q_norm * m_cos_theta);
        fout << omega << "\t" <<
            data->qz(ipt) << "\t" <<
            m_scale * data->reI(ipt) << "\t" <<
            m_scale * data->imI(ipt) << std::endl;
    }
    fout.close();

    std::cout << "Saved:\t" << data->m_nbSteps << std::endl;
}

bool RootEngine::done(const MCCalculator::Data * data)
{
    static bool isDone = false;

    if(data->m_imIMax != 0)
    {
        if(fabs(data->m_imIMax / data->m_reIMax) < m_programSettings.getEngineSettings().precision)
        {
            isDone = true;
        }
    }
    return isDone;
}


void RootEngine::copySettings(const std::string& orig_cfg, const std::string& out)
{
	boost::filesystem::path backup_cfg;

	backup_cfg = out;
	backup_cfg.replace_extension(".cfg");

	/* copy configuration 
	 * in a file with extension cfg and same filename as output file
	 */
	boost::filesystem::copy(orig_cfg, backup_cfg);
}

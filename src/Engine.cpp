#include "Engine.h"

Engine::Engine()
{
    m_sample = NULL;

    m_rng = NULL;
}

Engine::~Engine()
{
    if(m_sample)
    {
        delete m_sample;
    }
    if(m_rng)
    {
        gsl_rng_free (m_rng);
    }
}

void Engine::setup(const std::string& cfgfile, int seed)
{
    m_programSettings.read(cfgfile);

    setupRng(seed);
    setupSample();
}

MCCalculator * Engine::allocateMCCalculatorCoplanar()
{
    MCCalculator * calculator;
    double sigma_x, sigma_z;
    Geometry::Vector3d Q;
    MillerHexIndicesTransformator transformator(m_programSettings.getSampleSettings().a0,
                                                      m_programSettings.getSampleSettings().c0);

    Q = transformator.toVector3d(m_programSettings.getCalculatorSettings().Q);

    /**input resolutions are in q-units**/
    sigma_x = 1.0 / m_programSettings.getCalculatorSettings().qresolX;
    sigma_z = 1.0 / m_programSettings.getCalculatorSettings().qresolZ;

    calculator = new MCCalculatorCoplanar(m_rng, m_sample, Q, sigma_x, sigma_z);

    return calculator;
}

MCCalculator * Engine::allocateMCCalculatorSkewDouble()
{
    MCCalculator * calculator;
    double sigma_x, sigma_z;
    double lambda;
    Geometry::Vector3d Q;
    MillerHexIndicesTransformator transformator(m_programSettings.getSampleSettings().a0,
                                                      m_programSettings.getSampleSettings().c0);

    Q = transformator.toVector3d(m_programSettings.getCalculatorSettings().Q);
    /**input resolutions are in q-units**/
    sigma_x = 1.0 / m_programSettings.getCalculatorSettings().qresolX;
    sigma_z = 1.0 / m_programSettings.getCalculatorSettings().qresolZ;
    lambda = m_programSettings.getCalculatorSettings().lambda;

    calculator = new MCCalculatorSkewDouble(m_rng, m_sample, Q, sigma_x, sigma_z, lambda);

    return calculator;
}

MCCalculator * Engine::allocateMCCalculatorSkewTriple()
{
	MCCalculator * calculator;
    double sigma_x, sigma_z;
    double lambda;
    Geometry::Vector3d Q;
    MillerHexIndicesTransformator transformator(m_programSettings.getSampleSettings().a0,
                                                      m_programSettings.getSampleSettings().c0);

    Q = transformator.toVector3d(m_programSettings.getCalculatorSettings().Q);
    /**input resolutions are in q-units**/
    sigma_x = 1.0 / m_programSettings.getCalculatorSettings().qresolX;
    sigma_z = 1.0 / m_programSettings.getCalculatorSettings().qresolZ;
    lambda = m_programSettings.getCalculatorSettings().lambda;

    calculator = new MCCalculatorSkewTriple(m_rng, m_sample, Q, sigma_x, sigma_z, lambda);

    return calculator;
}

MCCalculator::Data * Engine::allocateMCCalculatorCoplanarData()
{
    MCCalculator::Data * data;
    size_t nbPoints, ipt;
    std::vector<double> qxvector, qzvector;

    nbPoints = m_programSettings.getEngineSettings().qxRange.m_sampling *
        m_programSettings.getEngineSettings().qzRange.m_sampling;

    data =
        new MCCalculator::Data(nbPoints);

    m_programSettings.getEngineSettings().qxRange.toVector(qxvector);
    m_programSettings.getEngineSettings().qzRange.toVector(qzvector);

    ipt = 0;
    for(size_t i = 0; i < qxvector.size(); ++i)
    {
        for(size_t j = 0; j < qzvector.size(); ++j)
        {
            data->qx(ipt) = qxvector[i];
            data->qz(ipt) = qzvector[j];
            data->reI(ipt) = 0.0;
            data->imI(ipt) = 0.0;
            ++ipt;
        }
    }
    data->m_nbSteps = m_programSettings.getEngineSettings().nbMCCalls;
    return data;
}

MCCalculator::Data * Engine::allocateMCCalculatorSkewTripleData()
{
    MCCalculator::Data * data;
    size_t nbPoints, ipt;
    std::vector<double> qxvector, qzvector;

    nbPoints = m_programSettings.getEngineSettings().qxRange.m_sampling *
        m_programSettings.getEngineSettings().qzRange.m_sampling;

    data =
        new MCCalculator::Data(nbPoints);

    m_programSettings.getEngineSettings().qxRange.toVector(qxvector);
    m_programSettings.getEngineSettings().qzRange.toVector(qzvector);

    ipt = 0;
    for(size_t i = 0; i < qxvector.size(); ++i)
    {
        for(size_t j = 0; j < qzvector.size(); ++j)
        {
            data->qx(ipt) = qxvector[i];
            data->qz(ipt) = qzvector[j];
            data->reI(ipt) = 0.0;
            data->imI(ipt) = 0.0;
            ++ipt;
        }
    }
    data->m_nbSteps = m_programSettings.getEngineSettings().nbMCCalls;
    return data;
}

MCCalculator::Data * Engine::allocateMCCalculatorSkewDoubleData()
{
    MCCalculator::Data * data;
    size_t nbPoints, ipt;
    std::vector<double> qzvector;

    nbPoints = m_programSettings.getEngineSettings().qzRange.m_sampling;

    data =
        new MCCalculator::Data(nbPoints);

    m_programSettings.getEngineSettings().qzRange.toVector(qzvector);

    ipt = 0;
    for(size_t j = 0; j < qzvector.size(); ++j)
    {
        data->qx(ipt) = 0.0;
        data->qz(ipt) = qzvector[j];
        data->reI(ipt) = 0.0;
        data->imI(ipt) = 0.0;
        ++ipt;
    }
    data->m_nbSteps = m_programSettings.getEngineSettings().nbMCCalls;
    return data;
}

void Engine::setupRng(int id)
{
    unsigned long int seed;

    m_rng = gsl_rng_alloc (gsl_rng_mt19937);
    seed = time(0) + id;
    gsl_rng_set (m_rng, seed);
}

void Engine::setupSample()
{
    MillerHexIndicesTransformator transformator(m_programSettings.getSampleSettings().a0,
                                                      m_programSettings.getSampleSettings().c0);
    /*double Qx, Qz;
    const MillerReciprocalHexIndices& Q = m_programSettings.getCalculatorSettings().Q;

    Qx = 2 * M_PI * sqrt(2.0 / 3 * (Q.H * Q.H + Q.K * Q.K + Q.I * Q.I))
	  / m_programSettings.getSampleSettings().a0;
	Qz = 2 * M_PI * Q.L / m_programSettings.getSampleSettings().c0;*/
    m_sample = new MCSampleHex(m_programSettings.getSampleSettings().thickness,
                                m_programSettings.getSampleSettings().width,
                               m_programSettings.getSampleSettings().nu,
                               Geometry::Vector3d(0.0, 0.0, 1.0),
                               m_programSettings.getSampleSettings().isHalfSpace);
    for(size_t i = 0; i < m_programSettings.getSampleSettings().misfit_interface.size(); ++i)
    {
        /**coefficient 1e-7 transforms density cm-1 -> nm-1**/
        m_sample->addMisfitInterface(m_programSettings.getSampleSettings().misfit_interface[i].rho * 1e-7,
                                     m_programSettings.getSampleSettings().misfit_interface[i].b_x,
                                     m_programSettings.getSampleSettings().misfit_interface[i].b_z,
                                     m_programSettings.getSampleSettings().thickness);
    }
    for(size_t ifam = 0; ifam < m_programSettings.getSampleSettings().threading_family.size(); ++ifam)
    {
        Geometry::Vector3d burgers;
        size_t nb_burg;
        double rc;
        double rho_fam, rho_burg;

        nb_burg = m_programSettings.getSampleSettings().threading_family[ifam].bs.size();
        rc = m_programSettings.getSampleSettings().threading_family[ifam].rc;
        rho_fam = m_programSettings.getSampleSettings().threading_family[ifam].rho;
        /**coefficient 1e-14 transforms density cm-2 -> nm-2**/
        rho_burg = rho_fam * 1e-14 / nb_burg;
        for(size_t iburg = 0; iburg < nb_burg; ++iburg)
        {
            burgers = transformator.toVector3d(m_programSettings.getSampleSettings().threading_family[ifam].bs[iburg]);

            m_sample->addThreadingLayer(m_rng,
                                    m_programSettings.getSampleSettings().thickness,
                                    rho_burg ,
                                    burgers,
                                    rc);
        }
    }
    for(size_t ifam = 0; ifam < m_programSettings.getSampleSettings().misfit_family.size(); ++ifam)
    {
        Geometry::Vector3d burgers, line;
        size_t nb_burg;
        int gamma;
        double rho_fam, rho_burg;

        nb_burg = m_programSettings.getSampleSettings().misfit_family[ifam].bs.size();
        gamma = m_programSettings.getSampleSettings().misfit_family[ifam].gamma;
        rho_fam = m_programSettings.getSampleSettings().misfit_family[ifam].rho;
        /**coefficient 1e-14 transforms density cm-2 -> nm-2**/
        rho_burg = rho_fam * 1e-7 / nb_burg;

        for(size_t iburg = 0; iburg < nb_burg; ++iburg)
        {
            burgers = transformator.toVector3d(m_programSettings.getSampleSettings().misfit_family[ifam].bs[iburg]);
            line = transformator.toVector3d(m_programSettings.getSampleSettings().misfit_family[ifam].ls[iburg]);

            m_sample->addMisfitSet(m_rng,
                                    m_programSettings.getSampleSettings().thickness,
                                    rho_burg,
                                    burgers,
                                    line,
                                    gamma);
        }
    }
}

void Engine::freeMCCalculator(MCCalculator * calculator)
{
    delete calculator;
}

void Engine::freeMCCalculatorData(MCCalculator::Data * data)
{
    delete data;
}

MCCalculator::Data * Engine::allocateMCCalculatorData()
{
    MCCalculator::Data * data;
    switch(m_programSettings.getEngineSettings().m_geometry)
    {
    case ProgramSettings::EngineSettings::geomCOPLANAR:
        data = allocateMCCalculatorCoplanarData();
        break;
    case ProgramSettings::EngineSettings::geomSKEW:
        if(m_programSettings.getEngineSettings().m_diffractometry == ProgramSettings::EngineSettings::diffDOUBLE)
        {
            data = allocateMCCalculatorSkewDoubleData();
        } else if(m_programSettings.getEngineSettings().m_diffractometry == ProgramSettings::EngineSettings::diffTRIPLE)
        {
			data = allocateMCCalculatorSkewTripleData();
        }
        else
        {
            data = NULL;
        }
        break;
    default:
        data = NULL;
        break;
    }
    return data;
}

MCCalculator * Engine::allocateMCCalculator()
{
    MCCalculator * calculator;
    switch(m_programSettings.getEngineSettings().m_geometry)
    {
    case ProgramSettings::EngineSettings::geomCOPLANAR:
        calculator = allocateMCCalculatorCoplanar();
        break;
    case ProgramSettings::EngineSettings::geomSKEW:
        if(m_programSettings.getEngineSettings().m_diffractometry == ProgramSettings::EngineSettings::diffDOUBLE)
        {
            calculator = allocateMCCalculatorSkewDouble();
        } else if(m_programSettings.getEngineSettings().m_diffractometry == ProgramSettings::EngineSettings::diffTRIPLE)
        {
			calculator = allocateMCCalculatorSkewTriple();
        }
        else
        {
            calculator = NULL;
        }
        break;
    default:
        calculator = NULL;
        break;
    }
    return calculator;
}

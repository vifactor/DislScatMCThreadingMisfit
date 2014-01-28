#include "MCCalculator.h"

MCCalculator::Data::Data(size_t nb_pts)
{
    m_nbPoints = nb_pts;
    m_nbSteps = 0;
    m_qx = new double[m_nbPoints];
    m_qz = new double[m_nbPoints];
    m_buffer = new double[2 * m_nbPoints];
    m_reIMax = 0;
    m_imIMax = 0;
}

MCCalculator::Data::~Data()
{
    delete[] m_qx;
    delete[] m_qz;
    delete[] m_buffer;
    m_nbPoints = 0;
    m_nbSteps = 0;
}

void MCCalculator::Data::append(const double * buff, unsigned long long nsteps)
{
    static double N0, N1, N;
    static double coef0, coef1;
    N0 = m_nbSteps;
    N1 = nsteps;
    N = N0 + N1;
    coef0 = N0 / N;
    coef1 = 1. / N;

    m_reIMax = coef0 * m_buffer[0] +
                        coef1 * buff[0];
    m_imIMax = coef0 * m_buffer[m_nbPoints] +
                        coef1 * buff[m_nbPoints];
    for(size_t ipt = 0; ipt < m_nbPoints; ++ipt)
    {
        m_buffer[ipt] = coef0 * m_buffer[ipt] +
                        coef1 * buff[ipt];
        m_buffer[ipt + m_nbPoints] = coef0 * m_buffer[ipt + m_nbPoints] +
                        coef1 * buff[ipt + m_nbPoints];

        if(m_buffer[ipt] > m_reIMax)
            m_reIMax = m_buffer[ipt];
        if(m_buffer[ipt + m_nbPoints] > m_imIMax)
            m_imIMax = m_buffer[ipt + m_nbPoints];
    }
    m_nbSteps = N;
}

void MCCalculator::Data::init()
{
    m_nbSteps = 0;
    for(size_t ipt = 0; ipt < 2 * m_nbPoints; ++ipt)
    {
        m_buffer[ipt] = 0.0;
    }
}

void MCCalculator::Data::transferTo(double * buff)
{
    for(size_t ipt = 0; ipt < 2 * m_nbPoints; ++ipt)
    {
        buff[ipt] = m_buffer[ipt];
        m_buffer[ipt] = 0.0;
    }
}

MCCalculator::MCCalculator()
{
    //ctor
}

MCCalculator::~MCCalculator()
{
    //dtor
}

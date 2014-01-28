#include "MisfitSet.h"
double MisfitSet::m_width = 0.0;

MisfitSet::MisfitSet(gsl_rng * rng, double depth, double rho, const Geometry::Vector3d& normal,
                           const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, int gamma)
{
    m_rng = rng;
    m_gamma = gamma;

    setDislocations(depth, rho, normal, Burgers, Line);
}

MisfitSet::~MisfitSet()
{
    for(size_t i = 0; i < m_dislocations.size(); ++i)
    {
        delete m_dislocations[i];
    }
}

const Geometry::Vector3d& MisfitSet::U(const Geometry::Vector3d& r) const
{
    m_u.set(0.0, 0.0, 0.0);
    for(size_t i = 0; i < m_dislocations.size(); ++i)
    {
        m_u += m_dislocations[i]->U(r);
    }

    return m_u;
}

void MisfitSet::update()
{
    static double pos;
    static double step;
    static int sx, sy, sz;

    /*sign of the misfit component does not change*/
    sx = 1;

    /*initial position is shifted by a mean distance to the left*/
    step = gsl_ran_gamma (m_rng, m_gamma, m_beta);
    pos = -m_width/2 - step;

    for(size_t i = 0; i < m_dislocations.size(); ++i)
    {
        /** setup new position */
        step = gsl_ran_gamma (m_rng, m_gamma, m_beta);
        pos += step;

        /** flip randomly by and bz components */
        sy = 1 - 2 * gsl_rng_uniform_int (m_rng, 2);
        sz = 1 - 2 * gsl_rng_uniform_int (m_rng, 2);

        m_dislocations[i]->moveTo(pos);
        m_dislocations[i]->multBurgers(sx, sy, sz);
    }
}

void MisfitSet::setDislocations(double depth, double rho, const Geometry::Vector3d& normal,
                           const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line)
{
    size_t nb;
    double mean;

    nb = rho * m_width;
    /*mean distance between dislocations*/
    if(nb > 1)
        mean = m_width / (nb - 1);
    else
        mean = m_width;

    m_beta = mean / m_gamma;

    m_dislocations.resize(nb);
    for(size_t i = 0; i < m_dislocations.size(); ++i)
    {
        m_dislocations[i] = new DislocationParallel(Burgers, Line, normal, depth);
    }
}

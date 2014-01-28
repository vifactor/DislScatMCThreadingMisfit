#include "ThreadingLayer.h"

double ThreadingLayer::m_width = 0.0;

ThreadingLayer::ThreadingLayer(gsl_rng * rng, double thickness, const Geometry::Vector3d& normal, double rho, const Geometry::Vector3d& Burgers, double Rc)
{
    m_rng = rng;
    m_thickness = thickness;
    m_Rc = Rc;
    setDislocations(normal, rho, Burgers);
}

ThreadingLayer::~ThreadingLayer()
{
    for(size_t i = 0; i < m_dislocations.size(); ++i)
    {
        delete m_dislocations[i];
    }
}

const Geometry::Vector3d& ThreadingLayer::U(const Geometry::Vector3d& r) const
{
    m_u.set(0.0, 0.0, 0.0);
    for(size_t i = 0; i < m_dislocations.size(); ++i)
    {
        m_u += m_dislocations[i]->U(r);
    }

    return m_u;
}

void ThreadingLayer::setDislocations(const Geometry::Vector3d& normal, double rho, const Geometry::Vector3d& Burgers)
{
    size_t nb;

    nb = rho * gsl_pow_2(m_width);

    if((nb % 2) != 0)
    {
        nb++;
    }
    m_dislocations.resize(nb);

    for(size_t i = 0; i < m_dislocations.size(); i += 2)
    {
        m_dislocations[i] = new DislocationPerpendicular(Burgers, normal);
        m_dislocations[i + 1] = new DislocationPerpendicular(-Burgers, normal);
    }
}

void ThreadingLayer::update()
{
    if(m_Rc > 0)
    {
        /** distribute threading dislocations as correlated pairs*/
        update_correlated();
    }
    else
    {
        /** distribute uncorrelated threading dislocations*/
        update_uncorrelated();
    }
}

void ThreadingLayer::update_correlated()
{
    static double R, phi, dx, dy;
    static Geometry::Vector2d center, pos1, pos2;

    for(size_t i = 0; i < m_dislocations.size(); i += 2)
    {
        center.x = gsl_ran_flat(m_rng, -m_width/2, m_width/2);
        center.y = gsl_ran_flat(m_rng, -m_width/2, m_width/2);

        R    = gsl_ran_gaussian_ziggurat(m_rng, m_Rc);
        phi  = gsl_ran_flat(m_rng, 0.0, 2.0 * M_PI);

        dx = R * cos(phi)/2;
        dy = R * sin(phi)/2;

        pos1.set(center.x - dx, center.y - dy);
        pos2.set(center.x + dx, center.y + dy);

        m_dislocations[i]->moveTo(pos1);
        m_dislocations[i + 1]->moveTo(pos2);
    }
}

void ThreadingLayer::update_uncorrelated()
{
    static Geometry::Vector2d pos;

    for(size_t i = 0; i < m_dislocations.size(); ++i)
    {
        pos.x = gsl_ran_flat(m_rng, -m_width/2, m_width/2);
        pos.y = gsl_ran_flat(m_rng, -m_width/2, m_width/2);

        m_dislocations[i]->moveTo(pos);
    }
}

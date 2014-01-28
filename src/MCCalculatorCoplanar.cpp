#include "MCCalculatorCoplanar.h"

MCCalculatorCoplanar::MCCalculatorCoplanar(gsl_rng * rng, MCSampleHex * sample, const Geometry::Vector3d& Q, double sigmax, double sigmaz)
{
    m_rng = rng;
    m_sample = sample;
    m_sigma_x = sigmax;
    m_sigma_z = sigmaz;
    m_Q = Q;

    setupLaboratoryFrame();
}

MCCalculatorCoplanar::~MCCalculatorCoplanar()
{
}

void MCCalculatorCoplanar::setupLaboratoryFrame()
{
    Geometry::Vector3d Qpar, Qper;

    Qpar = m_sample->m_normal * (m_sample->m_normal * m_Q);
    Qper = m_Q - Qpar;

    /*z along surface normal*/
    m_axis_z = m_sample->m_normal;
    m_axis_z.normalize();
    if(Qper.norm() <= m_epsilon)
    {
        /*for symmetric reflection orientation of x-axis is irrelevant*/
        /*find any vector perpendicular to z */
		Geometry::Vector3d vec;

		/*
		 * calculate the cross product [m_axis_z % vec],
		 * where vec is not collinear with m_axis_z
		 */
		if((m_axis_z.y != 0) || (m_axis_z.z != 0))
			m_axis_x = m_axis_z % Geometry::Vector3d(1, 0, 0);
		else
			m_axis_x = m_axis_z % Geometry::Vector3d(0, 1, 0);

		m_axis_x.normalize();
    }
    else
    {
        /*for asymmetric reflection orientation x along Qper*/
        m_axis_x = Qper.normalize();
    }
    /* y perpendicular to both */
    m_axis_y = m_axis_z % m_axis_x;

    /*std::cout << "Qx:\t" << m_Q * m_axis_x << std::endl;
    std::cout << "Qz:\t" << m_Q * m_axis_z << std::endl;*/


    /*initialize components of scattering vector in correlation function for misfit dislocations*/
    m_sample->setQ(m_Q * m_axis_x, m_Q * m_axis_y, m_Q * m_axis_z);

    std::cout << "Q_lab:\t" << m_Q * m_axis_x << "\t" << m_Q * m_axis_y << "\t" << m_Q * m_axis_z << std::endl;
}

void MCCalculatorCoplanar::add(Data * data)
{
    static double QU, qr, arg, Gmisfit;
    static double x, y, z1, z2, z;
    static Geometry::Vector3d r1, r2, U1, U2;

    /*generate two random points*/
    z1 = -gsl_ran_flat(m_rng, 0.0, m_sample->m_thickness);
    do
    {
        z = gsl_ran_gaussian_ziggurat(m_rng, m_sigma_z);
        z2 = z1 - z;
    } while((z2 < -m_sample->m_thickness) || (z2 > 0.0));

    x = gsl_ran_gaussian_ziggurat(m_rng, m_sigma_x);
    /*for coplanar geometry*/
    y = 0;

    /*find coordinates in the lab coord system*/
    /* minus sign before z-coordinates is due to the chosen laboratory coordinate frame*/
    r1 = m_axis_x * x / 2 + m_axis_z * z1;
    r2 = m_axis_x * (-x / 2) + m_axis_z * z2;
    U1 = m_sample->U_threading(r1);
    U2 = m_sample->U_threading(r2);

    QU = m_Q * (U1 - U2);
    Gmisfit = exp(-m_sample->T_misfit(-x, y, -z1, -z2));

    /*std::cout << "x, z1, z2, z:\t" <<  x << "\t" << z1 << "\t" << z2 << "\t" << z << std::endl;
    std::cout << "r1:\t" << r1.x << "\t" << r1.y << "\t" << r1.z << std::endl;
    std::cout << "r2:\t" << r2.x << "\t" << r2.y << "\t" << r2.z << std::endl;
    std::cout << "U1:\t" << U1.x << "\t" << U1.y << "\t" << U1.z << std::endl;
    std::cout << "U2:\t" << U2.x << "\t" << U2.y << "\t" << U2.z << std::endl;
    std::cout << "Q:\t" << m_Q.x << "\t" << m_Q.y << "\t" << m_Q.z << std::endl;
    std::cout << "QU:\t" << QU << std::endl;
    std::cout << "Gmisfit:\t" << Gmisfit << std::endl;*/

    for(size_t i = 0; i < data->size(); ++i)
    {
        qr = data->qx(i) * x + data->qz(i) * z;
        arg = QU + qr;

        data->reI(i) += cos(arg) * Gmisfit;
        data->imI(i) += sin(arg) * Gmisfit;
    }
}

void MCCalculatorCoplanar::run(Data * data)
{
    static size_t istep;

    istep = 0;
    while(istep < data->m_nbSteps)
    {
        m_sample->update();
        for(size_t idisl = 0; idisl < m_sample->getNbDislocations() + 1; ++idisl)
        {
            add(data);
            ++istep;
        }
    }
    data->m_nbSteps = istep;
}

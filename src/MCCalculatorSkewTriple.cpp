#include "MCCalculatorSkewTriple.h"

MCCalculatorSkewTriple::MCCalculatorSkewTriple(gsl_rng * rng, MCSampleHex * sample, const Geometry::Vector3d& Q, double sigmax, double sigmaz, double lambda)
{
    m_rng = rng;
    m_sample = sample;
    m_sigma_x = sigmax;
    m_sigma_z = sigmaz;
    m_Q = Q;
    m_lambda = lambda;

    setupLaboratoryFrame();
}

MCCalculatorSkewTriple::~MCCalculatorSkewTriple()
{
}

void MCCalculatorSkewTriple::setupLaboratoryFrame()
{
    Geometry::Vector3d Qpar, Qper;
    double cos_psi;

    /**laboratory coordinate frame for skew triple crystal is the same as for coplanar**/
    Qpar = m_sample->m_normal * inner_prod(m_sample->m_normal, m_Q);
    Qper = m_Q - Qpar;

    /*z along surface normal*/
    m_axis_z = m_sample->m_normal;
    normalize(m_axis_z);
    if(norm_2(Qper) <= m_epsilon)
    {
        /*for symmetric reflection orientation of x-axis is irrelevant*/
        /*find any vector perpendicular to z */
		Geometry::Vector3d vec;

		/*
		 * calculate the cross product [m_axis_z % vec],
		 * where vec is not collinear with m_axis_z
		 */
		if((m_axis_z[1] != 0) || (m_axis_z[2] != 0))
			m_axis_x = cross(m_axis_z, Geometry::Vector3d(1, 0, 0));
		else
			m_axis_x = cross(m_axis_z, Geometry::Vector3d(0, 1, 0));

		normalize(m_axis_x);
    }
    else
    {
        /*for asymmetric reflection orientation x along Qper*/
        m_axis_x = normalize(Qper);
    }
    /* y perpendicular to both */
    m_axis_y = cross(m_axis_z, m_axis_x);

    m_sin_psi = m_Q[2] / norm_2(m_Q);
    cos_psi = sqrt(1 - gsl_pow_2(m_sin_psi));

    m_cotan_psi = cos_psi / m_sin_psi;

    /*initialize components of scattering vector in correlation function for misfit dislocations*/
    m_sample->setQ(inner_prod(m_Q, m_axis_x),
                    inner_prod(m_Q, m_axis_y),
                    inner_prod(m_Q, m_axis_z));
}

void MCCalculatorSkewTriple::add(Data * data)
{
    static double QU, qr, arg, Gmisfit;
    static double x, y, z1, z2, z;
    static Geometry::Vector3d r1, r2, U1, U2;

    /*generate two random points*/
    z1 = -gsl_ran_flat(m_rng, 0.0, m_sample->m_thickness);
    if(fabs(m_sin_psi) > 0)
    {
        do
        {
            z = gsl_ran_gaussian_ziggurat(m_rng, m_sigma_z);
            z2 = z1 - z;
        } while((z2 < -m_sample->m_thickness) || (z2 > 0.0));
        x = z * m_cotan_psi;
        z = z / m_sin_psi;
    }
    else
    {
        z2 = z1;
        x = gsl_ran_gaussian_ziggurat(m_rng, m_sigma_x);
        z = x;
    }
    y = gsl_ran_gaussian_ziggurat(m_rng, m_sigma_x);

    /*find coordinates in the lab coord system*/
    /* minus sign before z-coordinates is due to the chosen laboratory coordinate frame*/
    r1 = m_axis_x * x / 2 + m_axis_y * y / 2 + m_axis_z * z1;
    r2 = m_axis_x * (-x / 2) + m_axis_y * (-y / 2) + m_axis_z * z2;
    U1 = m_sample->U_threading(r1);
    U2 = m_sample->U_threading(r2);

    QU = inner_prod(m_Q, (U1 - U2));
    /**
        the T correlation function is given in the following coordinate frame:
            - x-coordinate along in-plane Q
            - z1, z2- coordinates opposite to the surface normal
	**/

	/* !!!for the skew geometry with analyzer the correlation function should be T(x, y, z1, z2)!!! */
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
        qr = data->qx(i) * y + data->qz(i) * z;
        arg = QU + qr;

        data->reI(i) += cos(arg) * Gmisfit;
        data->imI(i) += sin(arg) * Gmisfit;
    }
}

void MCCalculatorSkewTriple::run(Data * data)
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


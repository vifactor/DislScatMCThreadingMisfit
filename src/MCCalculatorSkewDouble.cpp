#include "MCCalculatorSkewDouble.h"

using namespace boost;

MCCalculatorSkewDouble::MCCalculatorSkewDouble(gsl_rng * rng, MCSampleHex * sample, const Geometry::Vector3d& Q, double sigmax, double sigmaz, double lambda)
{
    m_rng = rng;
    m_sample = sample;
    m_sigma_x = sigmax;
    m_sigma_z = sigmaz;
    m_Q = Q;
    m_lambda = lambda;

    setupLaboratoryFrame();
}

MCCalculatorSkewDouble::~MCCalculatorSkewDouble()
{
}

void MCCalculatorSkewDouble::setupLaboratoryFrame()
{
    double theta;
    Geometry::Vector3d Q_cross_n, Kout, Kout_per, Kout_par, v1, v2;
    double cos_phi, sin_psi, cos_psi;

    theta = asin(norm_2(m_Q) * m_lambda/(4. * M_PI));

    v1 = m_Q / 2.0;
    /**vector perpendicular to normal and Q direction**/
    Q_cross_n = cross(m_Q, m_sample->m_normal);
    if(norm_2(Q_cross_n) == 0)
    {
        /*symmetric reflection: xy orientation is not important*/
        v2 = Geometry::Vector3d(1.0, 0.0, 0.0) * norm_2(v1) /tan(theta);
    }
    else
    {
        normalize(Q_cross_n);
        /** v2 = |Kout| cos(theta) Q_cross_n, |Kout = 2 pi / lambda| **/
        v2 = Q_cross_n * (2.0 * M_PI / m_lambda) * cos(theta);
    }
    Kout = v1 + v2;

    Kout_par = m_sample->m_normal * inner_prod(m_sample->m_normal, Kout);
    Kout_per = Kout - Kout_par;

    /*x along in-plane component of Kout*/
    m_axis_x = Kout_per;
    normalize(m_axis_x);

    /*z along surface normal*/
    m_axis_z = m_sample->m_normal;
    normalize(m_axis_z);

    /*z along surface normal*/
    m_axis_z = m_sample->m_normal;

    /*y perpendicular to both x and z*/
    m_axis_y = cross(m_axis_z, m_axis_x);

    //Phi  - angle between x-axis and Kout
    cos_phi = inner_prod(m_axis_x, Kout) / norm_2(Kout);
    m_sin_phi = sqrt(1.0 - gsl_pow_2(cos_phi));
    m_cotan_phi = cos_phi / m_sin_phi;

    //x axis along projection of Q on the surface
    sin_psi = m_Q[2] / norm_2(m_Q);
    /**
        limited precision problem required this type of hack
        otherwise for symmetric reflection when sin_psi = 0
        1 - sin_psi^2 appears to be a small negative double
        as a result cos_psi = nan (as well as cos_alpha = nan)
    */
    cos_psi = 1.0 - gsl_pow_2(sin_psi);
    if(cos_psi < 0)
    {
        cos_psi = 0.0;
    }
    else
    {
        cos_psi = sqrt(cos_psi);
    }

    m_cos_alpha = sin(theta) * cos_psi / cos_phi;

    /*std::cout << "sin_psi:\t" << sin_psi << std::endl;
    std::cout << "cos_psi:\t" << cos_psi << std::endl;
    std::cout << "cos_phi:\t" << cos_phi << std::endl;
    std::cout << "Phi:\t" << m_sin_phi << "\t" << m_cotan_phi << std::endl;
    std::cout << "alpha:\t" << m_cos_alpha << std::endl;

    std::cout << "axisX:\t" << m_axis_x.x << "\t" << m_axis_x.y << "\t" << m_axis_x.z << std::endl;
    std::cout << "axisY:\t" << m_axis_y.x << "\t" << m_axis_y.y << "\t" << m_axis_y.z << std::endl;
    std::cout << "axisZ:\t" << m_axis_z.x << "\t" << m_axis_z.y << "\t" << m_axis_z.z << std::endl;
    std::cout << m_Q.x << "\t" << m_Q.y << "\t" << m_Q.z << std::endl;*/

    /*initialize components of scattering vector in correlation function for misfit dislocations*/
    m_sample->setQ(inner_prod(m_Q, m_axis_x),
                    inner_prod(m_Q, m_axis_y),
                    inner_prod(m_Q, m_axis_z));

    std::cout << "Q_lab:\t" << inner_prod(m_Q, m_axis_x) 
            << "\t" << inner_prod(m_Q, m_axis_y) 
            << "\t" << inner_prod(m_Q, m_axis_z) << std::endl;
}

void MCCalculatorSkewDouble::add(Data * data)
{
    static double QU, qr, arg, Gmisfit;
    static double x, y, z1, z2, z;
    static Geometry::Vector3d r1, r2, U1, U2;

    /*generate two random points*/
    z1 = -gsl_ran_flat(m_rng, 0.0, m_sample->m_thickness);
    if(fabs(m_sin_phi) > 0)
    {
        do
        {
            z = gsl_ran_gaussian_ziggurat(m_rng, m_sigma_z);
            z2 = z1 - z;
        } while((z2 < -m_sample->m_thickness) || (z2 > 0.0));
        x = z * m_cotan_phi;
        z = z / m_sin_phi;
    }
    else
    {
        z2 = z1;
        x = gsl_ran_gaussian_ziggurat(m_rng, m_sigma_x);
        z = x;
    }
    /*for skew double geometry*/
    y = 0;

    /*find coordinates in the lab coord system*/
    /* minus sign before z-coordinates is due to the chosen laboratory coordinate frame*/
    r1 = m_axis_x * x / 2 + m_axis_z * z1;
    r2 = m_axis_x * (-x / 2) + m_axis_z * z2;
    U1 = m_sample->U_threading(r1);
    U2 = m_sample->U_threading(r2);

    QU = inner_prod(m_Q, (U1 - U2));
    /**
        the T correlation function is given in the following coordinate frame:
            - x-coordinate along in-plane Q
            - z1, z2- coordinates opposite to the surface normal**/
    Gmisfit = exp(-m_sample->T_misfit(-x * m_cos_alpha, y, -z1, -z2));

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
        qr = data->qz(i) * z;
        arg = QU + qr;

        data->reI(i) += cos(arg) * Gmisfit;
        data->imI(i) += sin(arg) * Gmisfit;
    }
}

void MCCalculatorSkewDouble::run(Data * data)
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


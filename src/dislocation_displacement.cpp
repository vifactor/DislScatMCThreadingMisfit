#include "dislocation_displacement.h"

void straight_dislocation_inf(double x, double y, double bx, double bz, double nu, double& u, double& v, double& w)
{
    const static double epsilon = 1e-6;
    double ux, vx, wx;
    double uz, vz, wz;
    double r2;

    r2 = x * x + y * y;

    ux = 0.0; vx = 0.0; wx = 0.0;
    if(fabs(bx) > epsilon)
    {
        ux = atan2(y,x) + x*y/(2 * (1.-nu) * r2);
        vx =-((1.-2*nu)/(4.*(1-nu))*log(r2)+x * x/(2.*(1.-nu)*r2));
        wx = 0.0;
    }

    uz = 0.0; vz = 0.0; wz = 0.0;
    if(fabs(bz) > epsilon)
    {
        uz = 0.0;
        vz = 0.0;
        wz = atan2(y,x);
    }

    u = (bx  * ux + bz * uz) / (2 * M_PI);
    v = (bx  * vx + bz * vz) / (2 * M_PI);
    w = (bx  * wx + bz * wz) / (2 * M_PI);
}

void straight_dislocation_perpendicular_hs(double x, double y, double z, double bx, double bz, double nu, double& u, double& v, double& w)
{
    const static double epsilon = 1e-6;
    static double ux, vx, wx;
    static double uz, vz, wz;
    static double r, r2, z2, y2, r_m_z2;

    r2 = x * x + y * y + z * z;
    r = sqrt(r2);
    r_m_z2 = gsl_pow_2(r - z);
    z2 = gsl_pow_2(z);
    y2 = gsl_pow_2(y);

    ux = 0.0; vx = 0.0; wx = 0.0;
    if(fabs(bx) > epsilon)
    {
        /**main terms correspond to Landau&Lifshitz, vol.VII, page 157 **/
        /**relaxation terms are from Lothe, in Indenbom&Lothe, page 367-368**/
        ux =//main term
            atan2(y,x) + x*y/(2 * (1.-nu) * (r2 - z2))
            //relaxation term
            +nu /(2 * (1.- nu)) * (2*x*y*z/(r * r_m_z2) + (1. - 2 * nu) * x*y/r_m_z2);
        vx =-((1.-2*nu)/(4.*(1. - nu)) * log(r2 - z2)+ x * x/(2.*(1.-nu)*(r2 - z2)))
            //relaxation term
            +nu/(2 * (1.0-nu)) * ((1. - 2 * nu) * log(r-z) - (3. - 2 * nu) * z / (r-z)
            +(3.-2*nu) * y2 / r_m_z2 - 2 * y2 / (r * (r-z)));
        wx =//relaxation term
            nu/(1.-nu) * y * (1./r + (1. - 2 * nu) * 1.0/(r-z));
    }

    uz = 0.0; vz = 0.0; wz = 0.0;
    if(fabs(bz) > epsilon)
    {
        uz = y/(r-z);
        vz = -x/(r-z);
        wz = atan2(y,x);
    }

    u = (bx  * ux + bz * uz) / (2 * M_PI);
    v = (bx  * vx + bz * vz) / (2 * M_PI);
    w = (bx  * wx + bz * wz) / (2 * M_PI);
}

void straight_dislocation_parallel_hs(double x, double z, double bx, double by, double bz, double nu, double d, double& u, double& v, double& w)
{
    const static double epsilon = 1e-6;
    static double u1xBx, u2xBx, u3xBx;
    static double u1zBx, u2zBx, u3zBx;
    static double u1xBz, u2xBz, u3xBz;
    static double u1zBz, u2zBz, u3zBz;
    static double u1yBy, u2yBy;
    static double alpha;
    static double cplus, cminus, cplus2;

    alpha = 1.0 / (2 * (1.0 - nu));
    u = 0;
    v = 0;
    w = 0;

    /**usefull variables often met in the expressions below**/
    cplus = x * x + (z + d) * (z + d);
    cminus = x * x + (z - d) * (z - d);
    cplus2 = gsl_pow_2(cplus);

    if(fabs(bx) > epsilon)
    {
        u1xBx = -(atan2(z-d, x) + alpha * x * (z - d)/cminus);
        u1zBx = ((1.0 - alpha)/2 * log(cminus)+alpha * x * x /cminus);

        u2xBx = (atan2(z+d, x) + alpha * x * (z + d)/cplus);
        u2zBx = -((1.0 - alpha)/2 * log(cplus)+alpha * x * x/cplus);

        u3xBx = 2.0 * d * ((1.0-alpha)* x /cplus-2 * alpha * x * z * (z+d)/cplus2);
        u3zBx = -2.0 * d * ((z+d)/cplus + alpha * z * ((z+d) * (z+d)-x * x)/cplus2);

        u += bx / (2 * M_PI)*(u1xBx + u2xBx + u3xBx);
        w += bx / (2 * M_PI)*(u1zBx + u2zBx + u3zBx);
    }

    if(fabs(by) > epsilon)
    {
        u1yBy =  atan2(x, z - d);
        u2yBy = -atan2(x, z + d);

        v += by/(2 * M_PI) * (u1yBy + u2yBy);
    }

    if(fabs(bz) > epsilon)
    {
        u1xBz = -((1.0 - alpha)/2 * log(cminus)+alpha * (z-d)* (z-d)/cminus);
        u1zBz = (atan2(x, z-d) + alpha * x * (z - d)/cminus);

        u2xBz = ((1.0 - alpha)/2*log(cplus)+alpha * (z+d) * (z+d)/cplus);
        u2zBz = -(atan2(x, z+d) + alpha *x * (z + d)/cplus);

        u3xBz = -2 * d *((1-alpha)*(z+d)/cplus+
                       alpha * (2 * x * x * z+d * cplus)/cplus2);
        u3zBz = -2 * d * ((1-alpha) * x/ cplus+2*alpha*x*z*(z+d)/cplus2);

        u += bz/(2 * M_PI)*(u1xBz + u2xBz + u3xBz);
        w += bz/(2 * M_PI)*(u1zBz + u2zBz + u3zBz);
    }
}

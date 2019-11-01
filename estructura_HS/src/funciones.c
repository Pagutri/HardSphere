#include "libraries.h"

double calculate_ck(double eta, double k)
{
    double alpha = 0.0;
    double delta = 0.0;
    double gamma = 0.0;
    double common_denom = 0.0;
    double ck = 0.0;
    double pi = 4.0 * atan(1.0);

    common_denom = pow(1.0 - eta, 4.0);

    alpha = -pow(1.0 + 2.0 * eta, 2.0) / common_denom;
    delta = 6.0 * eta * pow(1.0 + 0.5 * eta, 2.0) / common_denom;
    gamma = -0.5 * eta * alpha;

    ck = 4.0 * pi * (alpha * pow(k, 3.0) * (sin(k) - k * cos(k)) + \
         24.0 * gamma + delta * pow(k, 2.0) * (2.0 * k * sin(k) - \
         (pow(k, 2.0) - 2.0) * cos(k) - 2.0) + gamma * ((pow(k, 3.0) \
         - 24.0 * k) * sin(k) - (pow(k, 4.0) - 12.0 * pow(k, 2.0) + \
         24.0) * cos(k))) / pow(k, 6.0);

    return(ck);
}

double calculate_sk(double eta, double ck)
{
    double rho = 0.0;
    double sk = 0.0;
    double pi = 4.0 * atan(1.0);

    rho = 6.0 * eta / pi;
    sk = 1.0 / (1.0 - rho * ck);

    return(sk);
}
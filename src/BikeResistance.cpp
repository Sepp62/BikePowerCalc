/*
    Name:       BikeResistance.cpp
    Created:	26.03.2021 23:36:35
    Author:     Bernd-PC2\bernd2
*/

#include <stdio.h>
#include "BikeResistance.h"

#define PRINT_DEBUG

// https://mini.fandom.com/de/wiki/Fahrrad/Technik/Rollversuche
// Speed v on position  x
// v(x) = W*tan(arccos(exp(B*(x-K2))))
// with 
//  par[0] = W  = sqr( Cr*g/B )
//  par[1] = K2 = -1/B * ln(cos(arctan(v0/W)))
//  par[2] = B  = (rho*cw*A)/(2*m)
// Gradients: g[0] --> d/dW, g[1] --> d/dK2, g[2] --> d/dB
// https://www.ableitungsrechner.net/
// g[0] = d/dW  = exp(-B*(x-K2))* sqr(1-exp(2*B*(x-K2)))
// g[1] = d/dK2 = B*W*exp(-B*(x-K2))/sqr(1-exp(2*B*(x-K2)))
// g[2] = d/dB  = - (W*(x-K2)*exp(-B*(x-K2))/sqr(1-exp(2*B*(x-K2)))

/*
 * motion function for bike
 */
double motion_func(LevMarq::Vector& par, double x, void* fdata)
{
    // v(x) = W* tan(arccos(exp(B * (x - K2))))
    double W  = par[0];
    double K2 = par[1];
    double B  = par[2];
    return W * tan(acos(exp(B * (x - K2))));
}

/*
 * Gradient function for motion equation
 */
void gradient(LevMarq::Vector& g, LevMarq::Vector& par, double x, void* fdata)
{
    double W  = par[0];
    double K2 = par[1];
    double B  = par[2];
    g[0] = exp(-B*(x-K2))* sqrt(1-exp(2*B*(x-K2)));          // d/dW
    g[1] = B*W*exp(-B*(x-K2))/sqrt(1-exp(2*B*(x-K2)));       // d/dK2  
    g[2] = -W*(x-K2)*exp(-B*(x-K2))/sqrt(1-exp(2*B*(x-K2))); // d/dB
}

double BikeResistance::calcRho(double altitude, double temp)
{
    double rho = 1.293 * exp( -1.293 * m_g * altitude / 101300) * 273 / (1 * temp + 273);
    rho = round( rho * 1000 + .5) / 1000;
    return rho; 
}

void BikeResistance::SetSystemParams(double mass, double airTemp, double altitude)
{
    m_mass = mass;  // mass rider + mass bike + mass fictive (mf) --> mf = J/r^2 with J = moment of inertia and r = radius, mf is simplified mass wheel
    m_airTemp = airTemp;
    m_altitude = altitude;
}

bool BikeResistance::CalcResistance(double& cR, double& cwA, const measurePoints& points)
{
    double start_B = (m_rho * m_start_cwA) / (2 * m_mass);
    double start_W = sqrt(m_start_cr * m_g / start_B);
    double start_K2 = -1 / start_B * log(cos(atan(points[0].y / start_W))); // v0 = t_data[0].y = highest speed

    // W, K2, B
    Vector params{ start_W, start_K2, start_B }; // Initial values of parameters

    int nIterations;
    levmarq_init(&m_lmstat);
    // m_lmstat.verbose=1;
    m_lmstat.max_it      = 500;
    m_lmstat.init_lambda = 0.5;
    m_lmstat.up_factor   = 2;
    m_lmstat.down_factor = 2;
    m_lmstat.target_derr = 0.00001;

    m_rho = calcRho( m_altitude, m_airTemp );

    nIterations = levmarq( params, points, NULL, &motion_func, &gradient, NULL, &m_lmstat );

#ifdef PRINT_DEBUG
    printf("rho: %f\r\n", m_rho );
    printf("**************** End of calculation ***********************\n");
    printf("N iterations: %d\n", nIterations);
#endif 

    double W  = params[0];
    double K2 = params[1];
    double B  = params[2];

    cwA = B * 2 * m_mass/m_rho; // air resistance coefficent * surface
    cR  = W * W * B / m_g;      // roll resistance coefficent

#ifdef PRINT_DEBUG
    printf("cR: %f, cwA: %f\r\n", cR, cwA);
#endif

    return nIterations < m_lmstat.max_it;
}
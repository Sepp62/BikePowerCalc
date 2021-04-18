/*
    Name:       BikeResistance
    Created:	26.03.2021 23:36:35
    Author:     Bernd-PC2\bernd2
*/

#ifndef BIKERESISTANCE_H
#define BIKERESISTANCE_H

#include "levmarq.h"
#include <math.h>

// callbacks
double motion_func(LevMarq::Vector& par, double x, void* fdata);
void   gradient(LevMarq::Vector& g, LevMarq::Vector& par, double x, void* fdata);

class BikeResistance : public LevMarq
{
public:
    void SetSystemParams( double mass, double airTemp, double altitude ); // mass Rider + mass Bike + mass wheel
    bool CalcResistance( double & cR, double & cwA, const measurePoints & points );

protected:
    double m_mass      = 120.0; // kg
    double m_rho       = 1.202; // kg/m^3
    double m_airTemp   = 18.0;  // ° Celsius
    double m_altitude  = 0.0;
    double m_start_cr  = 0.006; // alt: 0.008, 0.51
    double m_start_cwA = 0.35;  // air resistance coefficent * surface
    const double m_g   = 9.81;  // m/s^2

    LMstat  m_lmstat;

    double calcRho( double altitude, double temp);
};
#endif  // BIKERESISTANCE_H

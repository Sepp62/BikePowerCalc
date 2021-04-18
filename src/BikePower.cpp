/*
    Name:       BikePower.cpp
    Created:    26.03.2021 23:36:35
    Author:     Bernd Wokoeck
*/

// #include <stdio.h>
#include "BikePower.h"

double BikePower::calcRho(double altitude, double temp)
{
    double rho = 1.293 * exp(-1.293 * m_g * altitude / 101300) * 273 / (1 * temp + 273);
    rho = round(rho * 1000 + .5) / 1000;
    return rho;
}

double BikePower::CurrentPower( double speedMeterPerSecond, double inclinationPercent, double altitude, double airTemp, unsigned long timestamp )
{
    // rho 
    double rho = m_lastRho;
    if (timestamp > m_nextRhoTime)
    {
        // save CPU and calc rho once per minute only
        rho = m_lastRho = calcRho(altitude, airTemp);
        m_nextRhoTime = timestamp + 60000L;
    }

    // inclination percent to degrees
    double inclinationRad = atan(inclinationPercent/100.0);

    // forces F
    double Froll  = m_mass * m_g * m_cR;    // weight force * roll coefficient, neglect inclination
    double FAir   = 0.5 * rho * m_cwA * speedMeterPerSecond * speedMeterPerSecond;
    double Fele   = m_mass * m_g * sin(inclinationRad);

    m_lastSpeed = speedMeterPerSecond;
    m_lastTime  = timestamp;

    // Power, P = F * v
    return (Froll + FAir + Fele) * speedMeterPerSecond;
}

double BikePower::Energy(double distance, double timeSeconds, double startAltitude, double endAltitude, double airTemp)
{
    double deltaH  = endAltitude - startAltitude;
    double rho     = calcRho(std::abs(deltaH) / 2.0, airTemp);
    double speedMs = distance/timeSeconds;
    double Wele  = m_mass * m_g * deltaH; // W = m * g * h
    double WAir  = 0.5 * rho * m_cwA * speedMs * speedMs * distance; // W = F*s
    double Wroll = m_mass * m_g * m_cR * distance; // W = F * s

    return Wele + WAir + Wroll;
}

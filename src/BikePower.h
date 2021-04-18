/*
    Name:       BikePower
    Created:    26.03.2021 23:36:35
    Author:     Bernd Wokoeck
*/

#ifndef BIKEPOWER_H
#define BIKEPOWER_H

#include <math.h>
#include <cmath>

class BikePower
{
protected:
    double m_mass    = 120.0;    // kg
    double m_cR      = 0.009725; // roll resistance - MTB on forest road https://www.leifiphysik.de/mechanik/reibung-und-fortbewegung/ausblick/reibungskraefte-beim-fahrradfahren
    double m_cwA     = 0.437392; // air resistance - cw * surface, MTB
    const double m_g = 9.81;     // m/s^2

    double        m_lastRho = 0.0;
    unsigned long m_nextRhoTime = 0L;
    double        calcRho(double altitude, double temp);

    double        m_lastSpeed = 0.0;
    unsigned long m_lastTime = 0L;

public:
    void   SetSystemParams(double mass, double cR, double cwA ) { m_mass = mass; m_cR = cR; m_cwA = cwA; } // mass Rider + mass Bike + mass wheel
    double CurrentPower(double speedMeterPerSecond, double inclinationPercent, double altitude, double airTemp, unsigned long timestamp );
    double Energy(double distance, double timeSeconds, double startAltitude, double endAltitude, double airTemp);
};

#endif // BIKEPOWER_H
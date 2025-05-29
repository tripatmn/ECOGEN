#include "EosPerfectGas.h"
#include <cmath>
#include <iostream>

//***********************************************************************

EosPerfectGas::EosPerfectGas(std::vector<std::string>& nameParameterEos, int& number) :
    Eos(number)
{
    nameParameterEos.push_back("gamma");
    nameParameterEos.push_back("R");
}

//***********************************************************************

EosPerfectGas::~EosPerfectGas() {}

//***********************************************************************

void EosPerfectGas::assignParametersEos(std::string name, std::vector<double> parametersEos)
{
    m_name = name;
    assert(parametersEos.size() == 2);
    m_gamma = parametersEos[0];
    m_R = parametersEos[1];
    m_cv = m_R / (m_gamma - 1);
}

//***********************************************************************

double EosPerfectGas::computeTemperature(const double& density, const double& pressure) const
{
    return pressure / (density * m_R);
}

//***********************************************************************

double EosPerfectGas::computeEnergy(const double& density, const double& pressure) const
{
    return pressure / ((m_gamma - 1) * density);
}

//***********************************************************************

double EosPerfectGas::computePressure(const double& density, const double& energy) const
{
    return (m_gamma - 1) * density * energy;
}

//***********************************************************************

double EosPerfectGas::computeDensity(const double& pressure, const double& temperature) const
{
    return pressure / (m_R * temperature);
}

//***********************************************************************

double EosPerfectGas::computeSoundSpeed(const double& density, const double& pressure) const
{
    return std::sqrt(m_gamma * pressure / density);
}

//***********************************************************************

double EosPerfectGas::computeEntropy(const double& temperature, const double& pressure) const
{
    return m_cv * std::log(temperature) - m_cv * std::log(pressure);
}

//***********************************************************************

double EosPerfectGas::computePressureIsentropic(const double& initialPressure, const double& initialDensity, const double& finalDensity) const
{
    return initialPressure * std::pow(finalDensity / initialDensity, m_gamma);
}

//***********************************************************************

double EosPerfectGas::computePressureHugoniot(const double& initialPressure, const double& initialDensity, const double& finalDensity) const
{
    return initialPressure * ((m_gamma + 1) * finalDensity - (m_gamma - 1) * initialDensity) /
           ((m_gamma + 1) * initialDensity - (m_gamma - 1) * finalDensity);
}

//***********************************************************************

double EosPerfectGas::computeDensityIsentropic(const double& initialPressure, const double& initialDensity, const double& finalPressure, double* drhodp) const
{
    double finalDensity = initialDensity * std::pow(finalPressure / initialPressure, 1.0 / m_gamma);
    if (drhodp) *drhodp = (1.0 / m_gamma) * (finalDensity / finalPressure);
    return finalDensity;
}

//***********************************************************************

double EosPerfectGas::computeDensityHugoniot(const double& initialPressure, const double& initialDensity, const double& finalPressure, double* drhodp) const
{
    double num = (m_gamma + 1) * finalPressure + (m_gamma - 1) * initialPressure;
    double denom = (m_gamma - 1) * finalPressure + (m_gamma + 1) * initialPressure;
    double finalDensity = initialDensity * num / std::max(denom, 1e-8);
    
    if (drhodp) 
        *drhodp = initialDensity * 4 * m_gamma * initialPressure / std::max(denom * denom, 1e-8);
    
    return finalDensity;
}


//***********************************************************************

double EosPerfectGas::computeDensityPfinal(const double& initialPressure, const double& initialDensity, const double& finalPressure, double* drhodp) const
{
    double finalDensity = initialDensity * (finalPressure / initialPressure);
    if (drhodp) *drhodp = initialDensity / initialPressure;
    return finalDensity;
}

//***********************************************************************

void EosPerfectGas::sendSpecialMixtureEos(double& gamPinfOverGamMinusOne, double& eRef, double& oneOverGamMinusOne, double& covolume) const
{
    gamPinfOverGamMinusOne = 0.0;
    eRef = 0.0;
    oneOverGamMinusOne = 1.0 / (m_gamma - 1);
    covolume = 0.0;
}

//***********************************************************************

void EosPerfectGas::verifyAndModifyPressure(double& pressure) const
{
    if (pressure < 1.e-8) pressure = 1.e-8;
}

//***********************************************************************


#include "EosPR.h"
#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>
#include <string>
#include <cfenv>
#include <cstdio>
#include <cstdlib>


#include <cfenv>

// Enable floating-point exceptions via a static initializer.
struct EnableFPExceptions {
    EnableFPExceptions() {
        // Enable traps for invalid operations, division by zero and overflow.
        // (The flags available depend on your compiler and system.)
        feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    }
};

// Create an instance to enable FP exceptions when the module loads.
static EnableFPExceptions enableFPEx;




// Although your compiler warns about #pragma STDC FENV_ACCESS, it is harmless
#pragma STDC FENV_ACCESS ON

// Macro to check for NaN and abort if found.
#define CHECK_NAN(x) do { \
    if (std::isnan(x)) { \
        std::fprintf(stderr, "NaN detected at %s:%d in variable %s\n", __FILE__, __LINE__, #x); \
        std::abort(); \
    } \
} while(0)

//***********************************************************************
// Constructor: Initialize PR-EOS parameters
EosPR::EosPR(std::vector<std::string>& nameParameterEos, int& number) : Eos(number)
{
    nameParameterEos.push_back("Tc");
    nameParameterEos.push_back("Pc");
    nameParameterEos.push_back("omega");
}

//***********************************************************************
// Destructor
EosPR::~EosPR(){}

//***********************************************************************
// Assign parameters from XML input
void EosPR::assignParametersEos(std::string name, std::vector<double> parametersEos)
{
    assert(parametersEos.size() == 3); // Expecting Tc, Pc, omega
    m_name = name;
    Tc    = parametersEos[0];
    Pc    = parametersEos[1];
    omega = parametersEos[2];

    calculateCoefficients();
}

//***********************************************************************
// Compute Peng-Robinson Coefficients (a, b, kappa)
void EosPR::calculateCoefficients()
{
    const double R = 8.314;
    a = 0.45724 * (R * R * Tc * Tc) / Pc;
    b = 0.07780 * (R * Tc) / Pc;
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
}

//***********************************************************************
// Compute temperature using PR-EOS
double EosPR::computeTemperature(const double& density, const double& pressure) const
{
    const double R = 8.314;
    double V = 1.0 / density;  // Molar volume

    if (fabs(V - b) < 1e-12) {
        throw std::runtime_error("computeTemperature: V - b is too small; state may be unphysical.");
    }

    // Residual function for temperature:
    auto f = [&](double T) -> double {
        double alpha = pow(1 + kappa * (1 - sqrt(T / Tc)), 2);
        double term1 = R * T / (V - b);
        double term2 = a * alpha / (V * (V + b) + b * (V - b));
        return term1 - term2 - pressure;
    };

    // Derivative of the residual function:
    auto df = [&](double T) -> double {
        double sqrtTr = sqrt(T / Tc);
        double dalpha_dT = - (kappa / Tc) * (1 / sqrtTr) * 2 * (1 + kappa * (1 - sqrtTr));
        double dterm1 = R / (V - b);
        double dterm2 = a * dalpha_dT / (V * (V + b) + b * (V - b));
        return dterm1 - dterm2;
    };

    double T_guess = (pressure * (V - b)) / R;
    const double tol = 1e-6;
    const int maxIter = 100;
    double T_new = T_guess;

    for (int iter = 0; iter < maxIter; ++iter) {
        double f_val = f(T_new);
        double df_val = df(T_new);
        if (fabs(df_val) < tol) {
            throw std::runtime_error("computeTemperature: Derivative too small during iteration.");
        }
        double T_old = T_new;
        T_new = T_old - f_val / df_val;
        if (fabs(T_new - T_old) < tol) {
            CHECK_NAN(T_new);
            return T_new;
        }
    }
    throw std::runtime_error("computeTemperature: Newton-Raphson did not converge.");
}

//***********************************************************************
// Compute specific heat capacity at constant volume (Cv) for PR-EOS
// Modified version that accepts temperature as an argument.
double EosPR::computeSpecificHeatCv(const double& density, const double& pressure, const double& temperature) const
{
    const double R = 8.314;
    double Tr = temperature / Tc;
    double alpha = pow(1 + kappa * (1 - sqrt(Tr)), 2);
    double aT = a * alpha;
    double denom = 1.0 - b * density;
    if (fabs(denom) < 1e-12) {
        denom = (denom > 0) ? 1e-12 : -1e-12;
    }
    double departureCv = - (aT / (b * R)) * (1.0 / denom);
    double Cv = (5.0 / 2.0) * R + departureCv;
    CHECK_NAN(Cv);
    return Cv;
}

//***********************************************************************
// Compute the adiabatic index (gamma) for PR-EOS
// Modified version to use the passed temperature.
double EosPR::computeGamma(const double& density, const double& pressure, const double& temperature) const
{
    const double R = 8.314;
    double cv = computeSpecificHeatCv(density, pressure, temperature);
    if (cv <= 0) {
        cv = 1e-12;
    }
    double gammaVal = 1.0 + (R / cv);
    CHECK_NAN(gammaVal);
    return gammaVal;
}

//***********************************************************************
// Compute pressure using PR-EOS equation
double EosPR::computePressure(const double& density, const double& temperature) const
{
    const double R = 8.314;
    double Tr = temperature / Tc;
    double alpha = pow(1 + kappa * (1 - sqrt(Tr)), 2);
    double aT = a * alpha;
    double pressure = (R * temperature / (density - b)) - (aT / (density * (density + b) + b * (density - b)));
    pressure = std::max(std::min(pressure, 1.0e7), 1.0e5);
    CHECK_NAN(pressure);
    return pressure;
}

//***********************************************************************
// Compute density from pressure and temperature using cubic EOS solver
double EosPR::computeDensity(const double& pressure, const double& temperature) const
{
    const double R = 8.314;
    double alpha = pow(1 + kappa * (1 - sqrt(temperature / Tc)), 2);
    double aT = a * alpha;
    double A = (aT * pressure) / (R * R * temperature * temperature);
    double B = (b * pressure) / (R * temperature);
    double Z = solveCubicEOS(A, B);
    double density = pressure / (R * temperature * Z);
    CHECK_NAN(density);
    return density;
}

//***********************************************************************
// Compute sound speed for PR-EOS
double EosPR::computeSoundSpeed(const double& density, const double& pressure) const
{
    const double R = 8.314;
    double T = computeTemperature(density, pressure);
    double alpha = pow(1 + kappa * (1 - sqrt(T / Tc)), 2);
    double aT = a * alpha;
    double dpdrho = std::max((R * T / (density - b)) - (aT / pow(density, 2)), 1e-8);
    double c = sqrt(dpdrho);
    CHECK_NAN(c);
    return c;
}

//***********************************************************************
// Compute internal energy for PR-EOS
double EosPR::computeEnergy(const double& density, const double& pressure) const
{
    const double R = 8.314;
    double T = computeTemperature(density, pressure);
    double alpha = pow(1 + kappa * (1 - sqrt(T / Tc)), 2);
    double aT = a * alpha;
    double logArgument = (density + (1 + sqrt(2)) * b) / (density + (1 - sqrt(2)) * b);
    if (logArgument <= 0) {
        throw std::runtime_error("Error in computeEnergy: Log argument is non-positive.");
    }
    double departureEnergy = - (aT / b) * (1.0 / (1.0 - b * density))
                             - (aT / (2.0 * sqrt(2.0) * b)) * log(logArgument);
    double energy = (computeSpecificHeatCv(density, pressure, T) * T) + departureEnergy;
    CHECK_NAN(energy);
    return energy;
}

//***********************************************************************
// Compute entropy using PR-EOS
double EosPR::computeEntropy(const double& temperature, const double& pressure) const
{
    double density = computeDensity(pressure, temperature);
    double maxDensity = 1.0 / b - 1.0;
    if (density > maxDensity) {
        density = maxDensity;
    }
    
    double gammaVal = computeGamma(density, pressure, temperature);
    double Cv = computeSpecificHeatCv(density, pressure, temperature);
    
    double logArg = std::pow(temperature, gammaVal) / std::max(std::pow(pressure, gammaVal - 1.), epsilonAlphaNull);
    if (logArg <= 0) {
        throw std::runtime_error("Error in computeEntropy: Log argument is non-positive.");
    }
    
    double entropy = Cv * log(logArg) + m_sRef;
    CHECK_NAN(entropy);
    return entropy;
}

//***********************************************************************
// Compute enthalpy using PR-EOS
double EosPR::computeEnthalpy(const double& density, const double& pressure) const
{
    const double R = 8.314;
    double T = computeTemperature(density, pressure);
    double Tr = T / Tc;
    double alpha = pow(1 + kappa * (1 - sqrt(Tr)), 2);
    double aT = a * alpha;
    double logArgument = (density + (1 + sqrt(2)) * b) / (density + (1 - sqrt(2)) * b);
    if (logArgument <= 0) {
        throw std::runtime_error("Error in computeEnthalpy: Log argument is non-positive.");
    }
    double enthalpy = (computeSpecificHeatCv(density, pressure, T) * T) + (R * T) - (aT / b) * log(logArgument);
    CHECK_NAN(enthalpy);
    return enthalpy;
}

//***********************************************************************
// Compute interface sound speed for PR-EOS
double EosPR::computeInterfaceSoundSpeed(const double& density, const double& interfacePressure, const double& pressure) const
{
    double gammaVal = computeGamma(density, pressure, computeTemperature(density, pressure));
    double c_interface = sqrt(((gammaVal - 1.) * interfacePressure + pressure) / std::max(density, epsilonAlphaNull));
    CHECK_NAN(c_interface);
    return c_interface;
}

//***********************************************************************
// Compute acoustic impedance for PR-EOS
double EosPR::computeAcousticImpedance(const double& density, const double& pressure) const
{
    double gammaVal = computeGamma(density, pressure, computeTemperature(density, pressure));
    double impedance = sqrt(density * gammaVal * pressure);
    CHECK_NAN(impedance);
    return impedance;
}

//***********************************************************************
// Compute density times interface sound speed square
double EosPR::computeDensityTimesInterfaceSoundSpeedSquare(const double& density, const double& interfacePressure, const double& pressure) const
{
    double gammaVal = computeGamma(density, pressure, computeTemperature(density, pressure));
    double value = (gammaVal - 1.) * interfacePressure + pressure;
    CHECK_NAN(value);
    return value;
}

//***********************************************************************
// Compute pressure along an isentropic path
double EosPR::computePressureIsentropic(const double& initialPressure, const double& initialDensity, const double& finalDensity) const
{
    double gammaVal = computeGamma(initialDensity, initialPressure, computeTemperature(initialDensity, initialPressure));
    double p = initialPressure * pow(finalDensity / std::max(initialDensity, epsilonAlphaNull), gammaVal);
    CHECK_NAN(p);
    return p;
}

//***********************************************************************
// Compute density along an isentropic path
double EosPR::computeDensityIsentropic(const double& initialPressure, const double& initialDensity, const double& finalPressure, double* drhodp) const
{
    double gammaVal = computeGamma(initialDensity, initialPressure, computeTemperature(initialDensity, initialPressure));
    double finalDensity = initialDensity * std::pow(finalPressure / std::max(initialPressure, epsilonAlphaNull), 1.0 / gammaVal);
    if (drhodp != nullptr) {
        *drhodp = finalDensity / std::max((gammaVal * finalPressure), epsilonAlphaNull);
    }
    CHECK_NAN(finalDensity);
    return finalDensity;
}

//***********************************************************************
// Compute density along the Hugoniot curve
double EosPR::computeDensityHugoniot(const double& initialPressure, const double& initialDensity, const double& finalPressure, double* drhodp) const
{
    double gammaVal = computeGamma(initialDensity, initialPressure, computeTemperature(initialDensity, initialPressure));
    double num = (gammaVal + 1.) * finalPressure + (gammaVal - 1.) * initialPressure;
    double denom = (gammaVal - 1.) * finalPressure + (gammaVal + 1.) * initialPressure;
    double finalDensity = initialDensity * num / std::max(denom, epsilonAlphaNull);
    if (drhodp != nullptr) {
        *drhodp = initialDensity * 4. * gammaVal * initialPressure / std::max((denom * denom), epsilonAlphaNull);
    }
    CHECK_NAN(finalDensity);
    return finalDensity;
}

//***********************************************************************
// Compute density given final pressure
double EosPR::computeDensityPfinal(const double& initialPressure, const double& initialDensity, const double& finalPressure, double* drhodp) const
{
    double gammaVal = computeGamma(initialDensity, initialPressure, computeTemperature(initialDensity, initialPressure));
    double num = gammaVal * finalPressure;
    double denom = num + initialPressure - finalPressure;
    double finalDensity = initialDensity * num / std::max(denom, epsilonAlphaNull);
    if (drhodp != nullptr) {
        *drhodp = initialDensity * gammaVal * initialPressure / std::max((denom * denom), epsilonAlphaNull);
    }
    CHECK_NAN(finalDensity);
    return finalDensity;
}

//***********************************************************************
// Compute pressure along the Hugoniot curve
double EosPR::computePressureHugoniot(const double& initialPressure, const double& initialDensity, const double& finalDensity) const
{
    double gammaVal = computeGamma(initialDensity, initialPressure, computeTemperature(initialDensity, initialPressure));
    double pressure = initialPressure * ((gammaVal + 1.) * finalDensity - (gammaVal - 1.) * initialDensity) /
           std::max(((gammaVal + 1.) * initialDensity - (gammaVal - 1.) * finalDensity), epsilonAlphaNull);
    CHECK_NAN(pressure);
    return pressure;
}

//***********************************************************************
// Solve cubic equation for Z in PR-EOS using Newton-Raphson method
double EosPR::solveCubicEOS(double A, double B) const
{
    const double tol = 1e-6;
    const int maxIter = 100;
    double Z = 1.0; // Initial guess
    
    for (int iter = 0; iter < maxIter; ++iter) {
        double f_val = Z * Z * Z - (1 - B) * Z * Z + (A - 3 * B * B - 2 * B) * Z - (A * B - B * B - B * B * B);
        double df_val = 3 * Z * Z - 2 * (1 - B) * Z + (A - 3 * B * B - 2 * B);
        if (fabs(df_val) < tol) {
            throw std::runtime_error("solveCubicEOS: Derivative is too close to zero, Newton-Raphson method fails.");
        }
        double Z_new = Z - f_val / df_val;
        if (fabs(Z_new - Z) < tol) {
            CHECK_NAN(Z_new);
            return Z_new;
        }
        Z = Z_new;
    }
    
    throw std::runtime_error("solveCubicEOS: Newton-Raphson did not converge.");
}

//***********************************************************************
// Compute isentropic enthalpy
double EosPR::computeEnthalpyIsentropic(const double& initialPressure, const double& initialDensity, 
    const double& finalPressure, double* dhdp) const
{
    const double R = 8.314;
    double finalDensity, drho;
    finalDensity = computeDensityIsentropic(initialPressure, initialDensity, finalPressure, &drho);
    double T_final = computeTemperature(finalDensity, finalPressure);
    double Tr = T_final / Tc;
    double alpha = pow(1 + kappa * (1 - sqrt(Tr)), 2);
    double aT = a * alpha;
    double A_final = (aT * finalPressure) / (R * R * T_final * T_final);
    double B_final = (b * finalPressure) / (R * T_final);
    double Z_final = solveCubicEOS(A_final, B_final);
    double initialT = computeTemperature(initialDensity, initialPressure);
    double enth = (computeSpecificHeatCv(finalDensity, finalPressure, T_final) * T_final) + (Z_final * R * T_final)
                  - (aT / b) * log((finalDensity + (1 + sqrt(2)) * b) / (finalDensity + (1 - sqrt(2)) * b));
    if (dhdp != nullptr) {
        *dhdp = (enth - (computeSpecificHeatCv(initialDensity, initialPressure, initialT) * initialT)) / (finalPressure - initialPressure);
    }
    CHECK_NAN(enth);
    return enth;
}

//***********************************************************************
// Compute density at saturation conditions
double EosPR::computeDensitySaturation(const double& pressure, const double& Tsat, const double& dTsatdP, double* drhodp) const
{
    double cv_value = computeSpecificHeatCv(pressure / Tsat, pressure, Tsat);
    if (drhodp != nullptr) {
        *drhodp = (computeGamma(pressure / Tsat, pressure, Tsat) - 1.) * cv_value * Tsat 
                  - pressure * (computeGamma(pressure / Tsat, pressure, Tsat) - 1.) * cv_value * dTsatdP;
    }
    double rho = pressure / std::max(((computeGamma(pressure / Tsat, pressure, Tsat) - 1.) * cv_value * Tsat), epsilonAlphaNull);
    CHECK_NAN(rho);
    return rho;
}

//***********************************************************************
// Compute density-energy product at saturation
double EosPR::computeDensityEnergySaturation(const double& pressure, const double& rho, const double& drhodp, double* drhoedp) const
{
    double gammaVal = computeGamma(rho, pressure, computeTemperature(rho, pressure));
    if (drhoedp != nullptr) {  
        *drhoedp = 1. / (gammaVal - 1.) + drhodp * m_eRef; 
    }
    double rhoe = pressure / (gammaVal - 1.) + rho * m_eRef;
    CHECK_NAN(rhoe);
    return rhoe;
}

//***********************************************************************
// Send special mixture EOS parameters for PR-EOS
void EosPR::sendSpecialMixtureEos(double& gamPinfOverGamMinusOne, double& eRef, double& oneOverGamMinusOne, double& covolume) const
{
    double referenceDensity = this->Pc / (8.314 * this->Tc);  // Ideal gas approximation
    double referencePressure = this->Pc;  
    double gammaVal = this->computeGamma(referenceDensity, referencePressure, computeTemperature(referenceDensity, referencePressure));
    gamPinfOverGamMinusOne = this->a / (this->b * (gammaVal - 1.));  
    eRef = this->m_eRef;
    oneOverGamMinusOne = 1. / (gammaVal - 1.);  
    covolume = this->b;
}

//***********************************************************************
// Compute function vfpfh
double EosPR::vfpfh(const double& pressure, const double& enthalpy) const
{
    double gammaVal = computeGamma(pressure / enthalpy, pressure, computeTemperature(pressure / enthalpy, pressure));
    double value = (gammaVal - 1.) * (enthalpy - m_eRef) / std::max((gammaVal * pressure), epsilonAlphaNull);
    CHECK_NAN(value);
    return value;
}

//***********************************************************************
// Compute derivative of density with respect to pressure
double EosPR::dvdhcp(const double& pressure) const
{
    double gammaVal = computeGamma(pressure / epsilonAlphaNull, pressure, computeTemperature(pressure / epsilonAlphaNull, pressure));
    double result = (gammaVal - 1.) / gammaVal / std::max(pressure, epsilonAlphaNull);
    CHECK_NAN(result);
    return result;
}

//***********************************************************************
// Compute derivative of density with respect to temperature at constant pressure
double EosPR::drhodpcT(const double& pressure, const double& temperature) const
{
    double gammaVal = computeGamma(pressure / temperature, pressure, temperature);
    double result = 1. / std::max(((gammaVal - 1.) * computeSpecificHeatCv(pressure / temperature, pressure, temperature) * temperature), epsilonAlphaNull);
    CHECK_NAN(result);
    return result;
}

//***********************************************************************
// Compute derivative of pressure with respect to enthalpy
double EosPR::dvdpch(const double& pressure, const double& enthalpy) const
{
    double gammaVal = computeGamma(pressure / enthalpy, pressure, computeTemperature(pressure / enthalpy, pressure));
    double result = (1. - gammaVal) / gammaVal * (enthalpy - m_eRef) / std::max((pressure * pressure), epsilonAlphaNull);
    CHECK_NAN(result);
    return result;
}

//***********************************************************************
// Verify pressure values to prevent numerical issues
void EosPR::verifyPressure(const double& pressure, const std::string& message) const
{
    if (pressure < 1.e-8)
        errors.push_back(Errors(message + " : too low pressure in EosPR"));
}

//***********************************************************************
// Adjust pressure to avoid numerical instability
void EosPR::verifyAndModifyPressure(double& pressure) const
{
    if (pressure < 1.e-8)
        pressure = 1.e-8;
}

//***********************************************************************
// Verify and Correct Maximum Density for PR-EOS
void EosPR::verifyAndCorrectDensityMax(double& density) const
{
    if (density > 1.0 / b - 1.0) {
        density = 1.0 / b - 1.0;
    }
}

//***********************************************************************
// Verify and Correct Maximum Density for PR-EOS (mass, alpha, density correction)
void EosPR::verifyAndCorrectDensityMax(const double& mass, double& alpha, double& density) const
{
    if (density > (1.0 / b) - 1.0) {
        density = (1.0 / b) - 1.0;
        alpha = mass / density;
    }
}


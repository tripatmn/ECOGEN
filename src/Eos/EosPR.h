#ifndef EOSPR_H
#define EOSPR_H

#include "Eos.h"
#include <vector>
#include <string>

//! \class EosPR
//! \brief Implements the Peng-Robinson Equation of State (PR-EOS)
class EosPR : public Eos {
public:
    EosPR(std::vector<std::string>& nameParameterEos, int& number);
    virtual ~EosPR();

    //! Assign parameters (Tc, Pc, omega) for different gases
    virtual void assignParametersEos(std::string name, std::vector<double> parametersEos);

    //! Compute temperature from density and pressure
    virtual double computeTemperature(const double& density, const double& pressure) const override;

    //! Compute pressure from density and energy
    virtual double computePressure(const double& density, const double& energy) const override;

    //! Compute density from pressure and temperature
    virtual double computeDensity(const double& pressure, const double& temperature) const override;

    //! Compute specific heat capacity at constant volume (Cv) for PR-EOS
    double computeSpecificHeatCv(const double& density, const double& pressure, const double& temperature) const;


    //! Compute the adiabatic index (gamma) for PR-EOS
    double computeGamma(const double& density, const double& pressure, const double& temperature) const;

    //! Compute pressure along an isentropic path
    virtual double computePressureIsentropic(const double& initialPressure, const double& initialDensity,
                                         const double& finalDensity) const override;

    //! Compute pressure along the Hugoniot curve
    virtual double computePressureHugoniot(const double& initialPressure, const double& initialDensity,
                                       const double& finalDensity) const override;

    //! Compute internal energy
    virtual double computeEnergy(const double& density, const double& pressure) const override;

    //! Compute enthalpy
    virtual double computeEnthalpy(const double& density, const double& pressure) const;

    //! Compute entropy
    virtual double computeEntropy(const double& temperature, const double& pressure) const override;

    //! Compute sound speed
    virtual double computeSoundSpeed(const double& density, const double& pressure) const override;

    //! Compute isentropic density
    virtual double computeDensityIsentropic(const double& initialPressure, const double& initialDensity, 
                                            const double& finalPressure, double* drhodp = nullptr) const override;

    //! Compute Hugoniot density
    virtual double computeDensityHugoniot(const double& initialPressure, const double& initialDensity, 
                                          const double& finalPressure, double* drhodp = nullptr) const override;

    //! Compute density given final pressure
    virtual double computeDensityPfinal(const double& initialPressure, const double& initialDensity, 
                                        const double& finalPressure, double* drhodp = nullptr) const override;

    //! Compute isentropic enthalpy
    virtual double computeEnthalpyIsentropic(const double& initialPressure, const double& initialDensity, 
                                             const double& finalPressure, double* dhdp = nullptr) const override;

    //! Compute density at saturation conditions
    virtual double computeDensitySaturation(const double& pressure, const double& Tsat, 
                                            const double& dTsatdP, double* drhodp = nullptr) const override;

    //! Compute density-energy product at saturation
    virtual double computeDensityEnergySaturation(const double& pressure, const double& rho, 
                                                  const double& drhodp, double* drhoedp = nullptr) const override;

    //! Solve cubic equation for Z-factor in PR-EOS
    double solveCubicEOS(double A, double B) const;

    //! Send special mixture EOS parameters
    virtual void sendSpecialMixtureEos(double& gamPinfOverGamMinusOne, double& eRef, 
                                       double& oneOverGamMinusOne, double& covolume) const override;
    //! Compute function vfpfh
    virtual double vfpfh(const double& pressure, const double& enthalpy) const override;
    
    //! Compute derivative of pressure with respect to enthalpy
    virtual double dvdpch(const double& pressure, const double& enthalpy) const override;
    
    //! Compute derivative of density with respect to pressure
    virtual double dvdhcp(const double& pressure) const override;
    
    //! Compute derivative of density with respect to temperature at constant pressure
    virtual double drhodpcT(const double& pressure, const double& temperature) const override;
    
    //! Verify pressure values to prevent numerical issues
    virtual void verifyPressure(const double& pressure, const std::string& message) const override;
    
    //! Adjust pressure to avoid numerical instability
    virtual void verifyAndModifyPressure(double& pressure) const override;
    //! Compute interface sound speed
    virtual double computeInterfaceSoundSpeed(const double& density, const double& interfacePressure, const double& pressure) const override;
    
    //! Compute acoustic impedance
    virtual double computeAcousticImpedance(const double& density, const double& pressure) const override;
    
    //! Compute density times interface sound speed square
    virtual double computeDensityTimesInterfaceSoundSpeedSquare(const double& density, const double& interfacePressure, const double& pressure) const override;
    
    virtual void verifyAndCorrectDensityMax(double& density) const override;
    virtual void verifyAndCorrectDensityMax(const double& mass, double& alpha, double& density) const override;
private:
    double Tc, Pc, omega; // Critical temperature, pressure, and acentric factor
    double a, b, kappa;    // PR-EOS parameters
    double m_eRef;  // Reference energy
    double m_sRef;  // Reference entropy

    void calculateCoefficients(); // Compute "a", "b", and "kappa" using Tc, Pc, omega
};

#endif // EOSPR_H


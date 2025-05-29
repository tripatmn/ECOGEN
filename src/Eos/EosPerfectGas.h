#ifndef EOSPERFECTGAS_H
#define EOSPERFECTGAS_H

#include "Eos.h"

//! \class EosPerfectGas
//! \brief Class describing a perfect gas equation of state
class EosPerfectGas : public Eos
{
public:
    EosPerfectGas(std::vector<std::string>& nameParameterEos, int& number);
    virtual ~EosPerfectGas();

    //! \brief Assign parameters for the EOS
    virtual void assignParametersEos(std::string name, std::vector<double> parametersEos) override;

    //! \brief Compute temperature
    virtual double computeTemperature(const double& density, const double& pressure) const override;

    //! \brief Compute internal energy
    virtual double computeEnergy(const double& density, const double& pressure) const override;

    //! \brief Compute pressure
    virtual double computePressure(const double& density, const double& energy) const override;

    //! \brief Compute density
    virtual double computeDensity(const double& pressure, const double& temperature) const override;

    //! \brief Compute sound speed
    virtual double computeSoundSpeed(const double& density, const double& pressure) const override;

    //! \brief Compute entropy
    virtual double computeEntropy(const double& temperature, const double& pressure) const override;

    //! \brief Compute isentropic pressure
    virtual double computePressureIsentropic(const double& initialPressure, const double& initialDensity, const double& finalDensity) const override;

    //! \brief Compute pressure along the Hugoniot curve
    virtual double computePressureHugoniot(const double& initialPressure, const double& initialDensity, const double& finalDensity) const override;

    //! \brief Compute density along an isentropic path
    virtual double computeDensityIsentropic(const double& initialPressure, const double& initialDensity, const double& finalPressure, double* drhodp = nullptr) const override;

    //! \brief Compute density along the Hugoniot curve
    virtual double computeDensityHugoniot(const double& initialPressure, const double& initialDensity, const double& finalPressure, double* drhodp = nullptr) const override;

    //! \brief Compute density during relaxation step
    virtual double computeDensityPfinal(const double& initialPressure, const double& initialDensity, const double& finalPressure, double* drhodp = nullptr) const override;

    //! \brief Send specific values for mixture EOS based on Perfect Gas
    virtual void sendSpecialMixtureEos(double& gamPinfOverGamMinusOne, double& eRef, double& oneOverGamMinusOne, double& covolume) const override;

    //! \brief Modify the pressure if its value is too low
    virtual void verifyAndModifyPressure(double& pressure) const override;

    //! \brief Do nothing for Perfect Gas
    virtual void verifyAndCorrectDensityMax(const double& /*mass*/, double& /*alpha*/, double& /*density*/) const override {};

    //! \brief Do nothing for Perfect Gas
    virtual void verifyAndCorrectDensityMax(double& /*density*/) const override {};

    //! \brief Get EOS type
    virtual TypeEOS getType() const override { return TypeEOS::PERFECTGAS; };

    //! \brief Get gamma
    virtual const double& getGamma() const override { return m_gamma; };

    //! \brief Get specific gas constant
    virtual const double& getR() const { return m_R; };

    //! \brief Get specific heat capacity at constant volume
    virtual const double& getCv() const override { return m_cv; };

private:
    double m_gamma;  //!< Adiabatic exponent of the fluid
    double m_R;      //!< Specific gas constant of the fluid
    double m_cv;     //!< Specific heat capacity at constant volume
};

#endif // EOSPERFECTGAS_H


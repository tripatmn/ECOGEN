#ifndef RELAXATIONPTMU_H
#define RELAXATIONPTMU_H

#include "Relaxation.h"
#include "../libTierces/LSODA.h"

//! \class RelaxationPTMu
//! \brief Pressure-Temperature-Chemical Potential relaxation / Phase change
class RelaxationPTMu : public Relaxation
{
public:
    //! \brief Relaxation constructor from an XML format reading
    //! \details Reading data from XML file under the following format:
    //!          ex: <dataPTMu liquid="SG_waterLiq.xml" vapor="IG_waterVap.xml" gas="IG_air.xml"/>
    //! \param element          XML element to read for source term
    //! \param fileName         string name of read XML file
    RelaxationPTMu(tinyxml2::XMLElement* element, std::vector<std::string> const& nameEOS, std::string fileName = "Unknown file");
    virtual ~RelaxationPTMu();

    //! \brief Initialize the theoretical critical pressure of the fluid
    //! \param cell           cell to get the eos
    //! \param numberPhases   number of phases
    virtual void initializeCriticalPressure(Cell *cell);

    //! \brief Stiff Thermo-Chemical relaxation method
    //! \details call for this method computes the thermodynamical equilibrium state in a given cell for a liquid, its vapor, and gas. Relaxed state is stored depending on the type enum
    //! \param cell           cell to relax
    //! \param dt             time step (not used here)
    //! \param type           enumeration allowing to relax either state in the cell or second order half time step state
    virtual void relaxation(Cell* cell, const double& /*dt*/, Prim type = vecPhases);

    //! \brief Return the pressure-, temperature- and chemical-potential-relaxation type
    virtual int getType() const { return PTMU; }

private:
    int m_liq;      //!< Liquid phase number for phase change
    int m_vap;      //!< Vapor phase number for phase change
    int m_gas;      //!< Gas phase number for phase interaction
    double m_pcrit; //!< Theoretical critical pressure of the fluid
};

#endif // RELAXATIONPTMU_H


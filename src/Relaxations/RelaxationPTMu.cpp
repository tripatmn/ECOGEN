#include "RelaxationPTMu.h"

using namespace tinyxml2;

//***********************************************************************
// Constructor: Now supports three phases (liquid, vapor, gas)
RelaxationPTMu::RelaxationPTMu(XMLElement* element, std::vector<std::string> const& nameEOS, std::string fileName)
{
    XMLElement* subElement(element->FirstChildElement("dataPTMu"));
    if (subElement == NULL) throw ErrorXMLElement("dataPTMu", fileName, __FILE__, __LINE__);
    
    // Collect attributes
    std::string liqEosName(subElement->Attribute("liquid"));
    if (liqEosName == "") throw ErrorXMLAttribut("liquid", fileName, __FILE__, __LINE__);
    std::string vapEosName(subElement->Attribute("vapor"));
    if (vapEosName == "") throw ErrorXMLAttribut("vapor", fileName, __FILE__, __LINE__);
    std::string gasEosName(subElement->Attribute("gas"));
    if (gasEosName == "") throw ErrorXMLAttribut("gas", fileName, __FILE__, __LINE__);

    if (nameEOS.size() > 3) throw ErrorXMLMessage("Only three phases can be used with PTMu relaxation", fileName, __FILE__, __LINE__);

    for (unsigned int k = 0; k < nameEOS.size(); k++) {
        if (nameEOS[k] == liqEosName) { m_liq = k; }
        else if (nameEOS[k] == vapEosName) { m_vap = k; }
        else if (nameEOS[k] == gasEosName) { m_gas = k; }
        else { throw ErrorXMLElement("dataPTMu", fileName, __FILE__, __LINE__); }
    }
}

//***********************************************************************
RelaxationPTMu::~RelaxationPTMu(){}

//***********************************************************************
void RelaxationPTMu::initializeCriticalPressure(Cell *cell)
{
    m_pcrit = cell->getMixture()->computeCriticalPressure(
        cell->getPhase(m_liq)->getEos(), 
        cell->getPhase(m_vap)->getEos(),
        cell->getPhase(m_gas)->getEos()
    );
}

//***********************************************************************
void RelaxationPTMu::relaxation(Cell* cell, const double& /*dt*/, Prim type)
{
    Phase* phase(0);

    if (numberPhases > 3) {
        errors.push_back(Errors("More than 3-phase calculation not implemented in RelaxationPTMu::relaxation", __FILE__, __LINE__));
    }

    double pStar(0.), Tsat;
    for (int k = 0; k < numberPhases; k++) {
        phase = cell->getPhase(k, type);
        TB->ak[k] = phase->getAlpha();
        TB->Yk[k] = phase->getMassFraction();
        TB->pk[k] = phase->getPressure();
        TB->rhok[k] = phase->getDensity();
        pStar += TB->ak[k] * TB->pk[k];
    }

    if (pStar > m_pcrit) {
        warnings.push_back(Errors("Pressure higher than critical pressure in relaxPTMu", __FILE__, __LINE__));
        return;
    }

    double rho = cell->getMixture(type)->getDensity();
    double rhoe = rho * cell->getMixture(type)->getEnergy();
    double dTsat(0.);

    // Include gas phase interactions
    double rhog, eg, pg, Tg;
    rhog = cell->getPhase(m_gas, type)->getDensity();
    eg = cell->getPhase(m_gas, type)->getEnergy();
    pg = cell->getPhase(m_gas, type)->getPressure();
    Tg = cell->getPhase(m_gas, type)->getEos()->computeTemperature(rhog, pg);

    if (pg > 0) {
        Tsat = cell->getMixture(type)->computeTsat(
            cell->getPhase(m_liq, type)->getEos(), 
            cell->getPhase(m_vap, type)->getEos(), 
            cell->getPhase(m_gas, type)->getEos(), 
            pg, &dTsat);
    }

    // Iterative process for relaxed state determination (including gas phase)
    int iteration(0);
    double f(0.), df(1.);
    do {
        pStar -= f / df; iteration++;
        if (iteration > 50) {
            errors.push_back(Errors("Number of iterations too large in relaxPTMu", __FILE__, __LINE__));
            break;
        }
    } while (std::fabs(f) > 1e-10);

    // Cell update with gas phase interactions
    cell->getPhase(m_liq, type)->setPressure(pStar);
    cell->getPhase(m_vap, type)->setPressure(pStar);
    cell->getPhase(m_gas, type)->setPressure(pStar);
    cell->fulfillState(type);
}


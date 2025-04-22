#include "Mixture.h"

int numberScalarsMixture;

//***************************************************************************

Mixture::Mixture(){}

//***************************************************************************

Mixture::~Mixture(){}

//***************************************************************************

void Mixture::printMixture(std::ofstream &fileStream) const
{
  //Scalar variables
  for (int var = 1; var <= this->getNumberScalars(); var++) {
    fileStream << this->returnScalar(var) << " ";
  }
  //Vector variables
  for (int var = 1; var <= this->getNumberVectors(); var++) {
    fileStream << this->returnVector(var).norm() << " ";
  }
}

//***************************************************************************

double Mixture::computeTsat(const Eos* eosLiq, const Eos* eosVap, const double& pressure, double* dTsat) {
  // Existing two-phase version
  if (!eosLiq || !eosVap) {
    Errors::errorMessage("Invalid EOS input to computeTsat");
    return 0.0;
  }

  double gammaL = eosLiq->getGamma();
  double pInfL = eosLiq->getPInf();
  double cvL = eosLiq->getCv();
  double e0L = eosLiq->getERef();
  double s0L = eosLiq->getSRef();

  double gammaV = eosVap->getGamma();
  double pInfV = eosVap->getPInf();
  double cvV = eosVap->getCv();
  double e0V = eosVap->getERef();
  double s0V = eosVap->getSRef();

  double A, B, C, D;
  A = (gammaL*cvL - gammaV*cvV + s0V - s0L) / (gammaV*cvV - cvV);
  B = (e0L - e0V) / (gammaV*cvV - cvV);
  C = (gammaV*cvV - gammaL*cvL) / (gammaV*cvV - cvV);
  D = (gammaL*cvL - cvL) / (gammaV*cvV - cvV);

  int iteration(0);
  double Tsat(0.1 * B / C);
  double f(0.), df(1.);
  do {
    Tsat -= f / df; iteration++;
    if (iteration > 50) {
      errors.push_back(Errors("Number of iterations too large in computeTsat", __FILE__, __LINE__));
      break;
    }
    f = A + B / Tsat + C * log(Tsat) - log(pressure + pInfV) + D * log(pressure + pInfL);
    df = C / Tsat - B / (Tsat * Tsat);
  } while (std::fabs(f) > 1e-10);

  double dfdp = -1. / (pressure + pInfV) + D / (pressure + pInfL);
  if (dTsat) *dTsat = -dfdp / df;

  return Tsat;
}

//***************************************************************************

double Mixture::computeTsat(const Eos* eosLiq, const Eos* eosVap, const Eos* eosGas, const double& pressure, double* dTsat) {
  // Three-phase version with gas pressure added to vapor term.
  if (!eosLiq || !eosVap || !eosGas) {
    Errors::errorMessage("Invalid EOS input to computeTsat");
    return 0.0;
  }

  double gammaL = eosLiq->getGamma();
  double pInfL = eosLiq->getPInf();
  double cvL = eosLiq->getCv();
  double e0L = eosLiq->getERef();
  double s0L = eosLiq->getSRef();

  double gammaV = eosVap->getGamma();
  double pInfV = eosVap->getPInf();
  double cvV = eosVap->getCv();
  double e0V = eosVap->getERef();
  double s0V = eosVap->getSRef();

  double gammaG = eosGas->getGamma();
  double pInfG = eosGas->getPInf();
  double cvG = eosGas->getCv();
  double e0G = eosGas->getERef();
  double s0G = eosGas->getSRef();

  double A, B, C, D;
  A = (gammaL * cvL - gammaV * cvV - gammaG * cvG + s0V - s0L + s0G) / (gammaV * cvV + gammaG * cvG - cvV - cvG);
  B = (e0L - e0V - e0G) / (gammaV * cvV + gammaG * cvG - cvV - cvG);
  C = (gammaV * cvV + gammaG * cvG - gammaL * cvL) / (gammaV * cvV + gammaG * cvG - cvV - cvG);
  // Modified D: now depends only on vapor properties.
  D = (gammaL * cvL - cvL) / (gammaV * cvV - cvV);

  int iteration(0);
  double Tsat(0.1 * B / C);
  double f(0.), df(1.);
  do {
    Tsat -= f / df;
    iteration++;
    if (iteration > 50) {
      errors.push_back(Errors("Number of iterations too large in computeTsat", __FILE__, __LINE__));
      break;
    }
    // Vapor term now uses (pressure + pInfV + pInfG), while liquid uses (pressure + pInfL)
    f = A + B / Tsat + C * log(Tsat) - log(pressure + pInfV + pInfG) + D * log(pressure + pInfL);
    df = C / Tsat - B / (Tsat * Tsat);
  } while (std::fabs(f) > 1e-10);

  double dfdp = -1. / (pressure + pInfV + pInfG) + D / (pressure + pInfL);
  if (dTsat) *dTsat = -dfdp / df;

  return Tsat;
}


//***************************************************************************

double Mixture::computePsat(const Eos* eosLiq, const Eos* eosVap, const double& temp)
{
  if (!eosLiq || !eosVap) {
    Errors::errorMessage("Invalid EOS input to computePsat");
    return 0.0;
  }

  double gammaL = eosLiq->getGamma();
  double pInfL = eosLiq->getPInf();
  double cvL = eosLiq->getCv();
  double e0L = eosLiq->getERef();
  double s0L = eosLiq->getSRef();

  double gammaV = eosVap->getGamma();
  double pInfV = eosVap->getPInf();
  double cvV = eosVap->getCv();
  double e0V = eosVap->getERef();
  double s0V = eosVap->getSRef();

  double A, B, C, D;
  A = (gammaL * cvL - gammaV * cvV + s0V - s0L) / (gammaV * cvV - cvV);
  B = (e0L - e0V) / (gammaV * cvV - cvV);
  C = (gammaV * cvV - gammaL * cvL) / (gammaV * cvV - cvV);
  D = (gammaL * cvL - cvL) / (gammaV * cvV - cvV);

  int iteration(0);
  double psat(2.e5);
  double f(0.), df(1.);
  do {
    psat -= f / df; iteration++;
    if (iteration > 50) {
      errors.push_back(Errors("Newton-Raphson did not converge in computePsat", __FILE__, __LINE__));
      break;
    }
    f = psat + pInfV - exp(A + B / temp + C * log(temp)) * std::pow(psat + pInfL, D);
    df = 1. - exp(A + B / temp + C * log(temp)) * D * std::pow(psat + pInfL, D - 1.);
  } while (std::fabs(f) > 1.e-8);

  return psat;
}

//***************************************************************************
// Three-phase version with gas pressure added to vapor term
double Mixture::computePsat(const Eos* eosLiq, const Eos* eosVap, const Eos* eosGas, const double& temp) {
  if (!eosLiq || !eosVap || !eosGas) {
    Errors::errorMessage("Invalid EOS input to computePsat");
    return 0.0;
  }

  double gammaL = eosLiq->getGamma();
  double pInfL = eosLiq->getPInf();
  double cvL = eosLiq->getCv();
  double e0L = eosLiq->getERef();
  double s0L = eosLiq->getSRef();

  double gammaV = eosVap->getGamma();
  double pInfV = eosVap->getPInf();
  double cvV = eosVap->getCv();
  double e0V = eosVap->getERef();
  double s0V = eosVap->getSRef();

  double gammaG = eosGas->getGamma();
  double pInfG = eosGas->getPInf();
  double cvG = eosGas->getCv();
  double e0G = eosGas->getERef();
  double s0G = eosGas->getSRef();

  double A, B, C, D;
  A = (gammaL * cvL - gammaV * cvV - gammaG * cvG + s0V - s0L + s0G) / (gammaV * cvV + gammaG * cvG - cvV - cvG);
  B = (e0L - e0V - e0G) / (gammaV * cvV + gammaG * cvG - cvV - cvG);
  C = (gammaV * cvV + gammaG * cvG - gammaL * cvL) / (gammaV * cvV + gammaG * cvG - cvV - cvG);
  // Modified D: remove gas term influence in the liquid component.
  D = (gammaL * cvL - cvL) / (gammaV * cvV - cvV);

  int iteration(0);
  double psat(2.e5);
  double f(0.), df(1.);
  do {
    psat -= f / df;
    iteration++;
    if (iteration > 50) {
      errors.push_back(Errors("Newton-Raphson did not converge in computePsat", __FILE__, __LINE__));
      break;
    }
    // Vapor term uses (psat + pInfV + pInfG), liquid remains (psat + pInfL)
    f = psat + pInfV + pInfG - exp(A + B / temp + C * log(temp)) * std::pow(psat + pInfL, D);
    df = 1. - exp(A + B / temp + C * log(temp)) * D * std::pow(psat + pInfL, D - 1.);
  } while (std::fabs(f) > 1.e-8);

  return psat;
}

//***************************************************************************

double Mixture::computeCriticalPressure(const Eos* eosLiq, const Eos* eosVap)
{
  double pCrit(0.);

  double gammaL = eosLiq->getGamma();
  double pInfL = eosLiq->getPInf();
  double cvL = eosLiq->getCv();

  double gammaV = eosVap->getGamma();
  double pInfV = eosVap->getPInf();
  double cvV = eosVap->getCv();

  // vVsat = vLsat at critical point
  pCrit = pInfL * (gammaV - 1.) * cvV - pInfV * (gammaL - 1.) * cvL;
  pCrit /= (gammaL - 1.) * cvL - (gammaV - 1.) * cvV;

  return pCrit;
}

//***************************************************************************

double Mixture::computeCriticalPressure(const Eos* eosLiq, const Eos* eosVap, const Eos* eosGas)
{
  
  double gammaL = eosLiq->getGamma();
  double pInfL  = eosLiq->getPInf();
  double cvL    = eosLiq->getCv();
  
  double gammaV = eosVap->getGamma();
  double pInfV  = eosVap->getPInf();
  double cvV    = eosVap->getCv();
  
  double pInfG  = eosGas->getPInf();
  
  // Effective vapor pressure offset: gas pressure is added to the vapor term.
  double pInf_eff = pInfV + pInfG;
  
  // Compute critical pressure analogous to the two-phase version using effective vapor properties.
  double pCrit = pInfL * (gammaV - 1.) * cvV - pInf_eff * (gammaL - 1.) * cvL;
  pCrit /= (gammaL - 1.) * cvL - (gammaV - 1.) * cvV;
  
  return pCrit;
}


//***************************************************************************


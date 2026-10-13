#include "GGEMS/materials/GGEMSMaterials.hh"
#include "GGEMS/physics/GGEMSPairProduction.hh"
#include "GGEMS/physics/GGEMSProcessConstants.hh"
#include "GGEMS/tools/GGEMSTools.hh"

GGEMSPairProduction::GGEMSPairProduction(std::string const& primary_particle, bool const& is_secondary)
: GGEMSEMProcess()
{
  GGcout("GGEMSPairProduction", "GGEMSPairProduction", 3) << "GGEMSPairProduction creating..." << GGendl;

  if (primary_particle != "gamma") {
    std::ostringstream oss(std::ostringstream::out);
    oss << "For Pair Production, incident particle has to be a 'gamma'";
    GGEMSMisc::ThrowException("GGEMSPairProduction", "GGEMSPairProduction", oss.str());
  }

  process_name_ = "PairProduction";
  process_id_ = PAIR_PRODUCTION;
  primary_particle_ = "gamma";
  secondary_particle_ = "e-";
  is_secondaries_ = is_secondary;

  GGcout("GGEMSPairProduction", "GGEMSPairProduction", 3) << "GGEMSPairProduction created!!!" << GGendl;
}

GGEMSPairProduction::~GGEMSPairProduction(void)
{
  GGcout("GGEMSPairProduction", "~GGEMSPairProduction", 3) << "GGEMSPairProduction erasing..." << GGendl;
  GGcout("GGEMSPairProduction", "~GGEMSPairProduction", 3) << "GGEMSPairProduction erased!!!" << GGendl;
}

GGfloat GGEMSPairProduction::ComputeCrossSectionPerAtom(GGfloat const& energy, GGuchar const& atomic_number) const
{
  static GGfloat const kElectronMassC2 = 0.510998910f;
  static GGfloat const kThreshold = 2.0f * kElectronMassC2;
  static GGfloat const kEnergyLimit = 1.5f;
  static GGfloat const kMicrobarn = 1.0e-28f;

  if (energy < kThreshold || atomic_number < 1) return 0.0f;

  GGfloat const Z = static_cast<GGfloat>(atomic_number);

  static GGfloat const a0 =  8.7842e+2f, a1 = -1.9625e+3f, a2 =  1.2949e+3f, a3 = -2.0028e+2f, a4 =  1.2575e+1f, a5 = -2.8333e-1f;
  static GGfloat const b0 = -1.0342e+1f, b1 =  1.7692e+1f, b2 = -8.2381f,    b3 =  1.3063f,    b4 = -9.0815e-2f, b5 =  2.3586e-3f;
  static GGfloat const c0 = -4.5263e+2f, c1 =  1.1161e+3f, c2 = -8.6749e+2f, c3 =  2.1773e+2f, c4 = -2.0467e+1f, c5 =  6.5372e-1f;

  GGfloat e_eff = energy < kEnergyLimit ? kEnergyLimit : energy;
  GGfloat X  = logf(e_eff / kElectronMassC2);
  GGfloat X2 = X*X, X3 = X2*X, X4 = X3*X, X5 = X4*X;

  GGfloat F1 = a0 + a1*X + a2*X2 + a3*X3 + a4*X4 + a5*X5;
  GGfloat F2 = b0 + b1*X + b2*X2 + b3*X3 + b4*X4 + b5*X5;
  GGfloat F3 = c0 + c1*X + c2*X2 + c3*X3 + c4*X4 + c5*X5;

  GGfloat xsection_microbarn = (Z + 1.0f) * (F1*Z + F2*Z*Z + F3);
  if (xsection_microbarn < 0.0f) xsection_microbarn = 0.0f;

  // below 1.5 MeV: scale to zero at threshold, as in Geant4 G4BetheHeitlerModel
  if (energy < kEnergyLimit) {
    GGfloat x = (energy - kThreshold) / (kEnergyLimit - kThreshold);
    xsection_microbarn *= x * x;
  }
  return xsection_microbarn * kMicrobarn;
}

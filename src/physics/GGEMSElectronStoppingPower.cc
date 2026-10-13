// ************************************************************************
// * This file is part of GGEMS.                                          *
// * Extension: Electron stopping-power table (NIST ESTAR based)          *
// * Added as part of thesis work extending GGEMS to electron transport   *
// ************************************************************************

/*!
  \file GGEMSElectronStoppingPower.cc

  \brief GGEMS class handling electron stopping-power tables (from NIST ESTAR data)
*/

#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <vector>

#include "GGEMS/physics/GGEMSElectronStoppingPower.hh"
#include "GGEMS/physics/GGEMSParticleCrossSections.hh"
#include "GGEMS/materials/GGEMSMaterialTables.hh"
#include "GGEMS/global/GGEMSOpenCLManager.hh"
#include "GGEMS/tools/GGEMSPrint.hh"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

GGEMSElectronStoppingPower::GGEMSElectronStoppingPower(void)
{
  GGcout("GGEMSElectronStoppingPower", "GGEMSElectronStoppingPower", 3) << "Allocation of GGEMSElectronStoppingPower..." << GGendl;
  process_name_ = "ElectronIonisation";
  primary_particle_ = "e-";
  secondary_particle_ = "";
  is_secondaries_ = false;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

GGEMSElectronStoppingPower::~GGEMSElectronStoppingPower(void)
{
  GGcout("GGEMSElectronStoppingPower", "~GGEMSElectronStoppingPower", 3) << "Deallocation of GGEMSElectronStoppingPower..." << GGendl;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void GGEMSElectronStoppingPower::LoadStoppingPowerTable(std::string const& material_name, std::string const& file_path)
{
  GGcout("GGEMSElectronStoppingPower", "LoadStoppingPowerTable", 1) << "Loading NIST ESTAR table for material '" << material_name << "' from " << file_path << GGendl;

  loaded_material_name_ = material_name;
  nist_energy_MeV_.clear();
  nist_stopping_power_.clear();

  std::ifstream infile(file_path);
  if (!infile.is_open()) {
    GGEMSMisc::ThrowException("GGEMSElectronStoppingPower", "LoadStoppingPowerTable", "Could not open stopping power file: " + file_path);
  }

  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty() || line[0] == '#') continue; // skip comments/blank lines

    std::istringstream iss(line);
    GGfloat energy, s_el, s_rad, s_tot, csda_range, rad_yield, delta;
    if (!(iss >> energy >> s_el >> s_rad >> s_tot >> csda_range >> rad_yield >> delta)) continue;

    nist_energy_MeV_.push_back(energy);
    nist_stopping_power_.push_back(s_tot);
  }
  infile.close();

  GGcout("GGEMSElectronStoppingPower", "LoadStoppingPowerTable", 1) << "Loaded " << nist_energy_MeV_.size() << " energy points" << GGendl;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

GGfloat GGEMSElectronStoppingPower::InterpolateStoppingPower(GGfloat const& energy_MeV) const
{
  if (nist_energy_MeV_.empty()) return 0.0f;

  // Clamp to table range
  if (energy_MeV <= nist_energy_MeV_.front()) return nist_stopping_power_.front();
  if (energy_MeV >= nist_energy_MeV_.back()) return nist_stopping_power_.back();

  // Find bracketing points
  auto it = std::lower_bound(nist_energy_MeV_.begin(), nist_energy_MeV_.end(), energy_MeV);
  size_t idx_high = std::distance(nist_energy_MeV_.begin(), it);
  size_t idx_low = idx_high - 1;

  GGfloat e_low = nist_energy_MeV_[idx_low];
  GGfloat e_high = nist_energy_MeV_[idx_high];
  GGfloat s_low = nist_stopping_power_[idx_low];
  GGfloat s_high = nist_stopping_power_[idx_high];

  // Log-log interpolation (stopping power vs energy is smoother in log-log space)
  GGfloat log_e = logf(energy_MeV);
  GGfloat log_e_low = logf(e_low);
  GGfloat log_e_high = logf(e_high);
  GGfloat log_s_low = logf(s_low);
  GGfloat log_s_high = logf(s_high);

  GGfloat frac = (log_e - log_e_low) / (log_e_high - log_e_low);
  GGfloat log_s = log_s_low + frac * (log_s_high - log_s_low);

  return expf(log_s);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

GGfloat GGEMSElectronStoppingPower::ComputeCrossSectionPerAtom(GGfloat const& /*energy*/, GGuchar const& /*atomic_number*/) const
{
  // Not used: this process is table-based (NIST ESTAR), not formula-based.
  // Required override from GGEMSEMProcess base class.
  return 0.0f;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

GGfloat GGEMSElectronStoppingPower::BetheStoppingPowerPerElement(GGfloat const& energy_MeV, GGint const& Z, GGfloat const& A, GGfloat const& I_eV) const
{
  // Bethe collision stopping power (ICRU-37 / Berger-Seltzer form), electrons
  const GGfloat mc2 = 0.510999f; // electron rest mass energy, MeV
  const GGfloat K = 0.153536f;   // MeV.cm2/g (2*pi*N_A*r_e^2*mc2, standard constant)

  GGfloat tau = energy_MeV / mc2;               // kinetic energy in units of mc2
  GGfloat beta2 = 1.0f - 1.0f / ((1.0f + tau) * (1.0f + tau));
  if (beta2 <= 0.0f) beta2 = 1.0e-6f; // guard against zero/negative at extremely low energy

  GGfloat I_MeV = I_eV * 1.0e-6f; // mean excitation energy, MeV

  // F-(tau) term for electrons (ICRU-37 eq. for the Moller/collision stopping power)
  GGfloat F_minus = 1.0f - beta2 + (tau * tau / 8.0f - (2.0f * tau + 1.0f) * logf(2.0f)) / ((tau + 1.0f) * (tau + 1.0f));

  GGfloat log_term = logf((tau * tau * (tau + 2.0f)) / 2.0f * (mc2 / I_MeV) * (mc2 / I_MeV));

  // Mass collision stopping power, MeV.cm2/g (no density-effect correction, v1 -- valid well below the
  // density-effect-significant energy range for these purposes)
  GGfloat stopping_power = (K / beta2) * (static_cast<GGfloat>(Z) / A) * (log_term + F_minus);

  return stopping_power;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace {
  // Elemental data lookup by atomic number Z: returns true and fills A, I if Z is known.
  // Covers elements present in GGEMS materials.txt / materials_icrp110.txt (biological tissues + common structural elements).
  bool LookupElementByZ(GGuchar const& Z, GGfloat& A, GGfloat& I_eV)
  {
    switch (Z) {
      case 1:  A = 1.008f;  I_eV = 19.2f;  return true; // H
      case 6:  A = 12.011f; I_eV = 78.0f;  return true; // C
      case 7:  A = 14.007f; I_eV = 82.0f;  return true; // N
      case 8:  A = 15.999f; I_eV = 95.0f;  return true; // O
      case 9:  A = 18.998f; I_eV = 115.0f; return true; // F
      case 11: A = 22.990f; I_eV = 149.0f; return true; // Na
      case 12: A = 24.305f; I_eV = 156.0f; return true; // Mg
      case 13: A = 26.982f; I_eV = 166.0f; return true; // Al
      case 14: A = 28.085f; I_eV = 173.0f; return true; // Si
      case 15: A = 30.974f; I_eV = 173.0f; return true; // P
      case 16: A = 32.06f;  I_eV = 180.0f; return true; // S
      case 17: A = 35.45f;  I_eV = 174.0f; return true; // Cl
      case 19: A = 39.098f; I_eV = 190.0f; return true; // K
      case 20: A = 40.078f; I_eV = 191.0f; return true; // Ca
      case 26: A = 55.845f; I_eV = 286.0f; return true; // Fe
      case 30: A = 65.38f;  I_eV = 330.0f; return true; // Zn
      case 53: A = 126.904f;I_eV = 491.0f; return true; // I
      default: return false;
    }
  }
}

GGfloat GGEMSElectronStoppingPower::ComputeBetheStoppingPower(GGfloat const& energy_MeV, std::vector<ElementComposition> const& composition) const
{
  // Bragg additivity: total mixture stopping power = weighted sum of elemental stopping powers
  GGfloat total_stopping_power = 0.0f;
  for (auto const& elem : composition) {
    GGfloat A, I_eV;
    if (!LookupElementByZ(elem.Z, A, I_eV)) {
      GGwarn("GGEMSElectronStoppingPower", "ComputeBetheStoppingPower", 0) << "Element Z=" << static_cast<int>(elem.Z) << " not in table, skipping (contributes 0)" << GGendl;
      continue;
    }
    GGfloat elem_sp = BetheStoppingPowerPerElement(energy_MeV, elem.Z, A, I_eV);
    total_stopping_power += elem.weight_fraction * elem_sp;
  }

  return total_stopping_power;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void GGEMSElectronStoppingPower::BuildCrossSectionTables(cl::Buffer* particle_cross_sections, cl::Buffer* material_tables, GGsize const& thread_index)
{
  GGcout("GGEMSElectronStoppingPower", "BuildCrossSectionTables", 1) << "Building electron stopping power table (multi-material, Bethe formula)..." << GGendl;

  GGEMSOpenCLManager& opencl_manager = GGEMSOpenCLManager::GetInstance();

  GGEMSParticleCrossSections* cs = opencl_manager.GetDeviceBuffer<GGEMSParticleCrossSections>(particle_cross_sections, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, sizeof(GGEMSParticleCrossSections), thread_index);
  GGEMSMaterialTables* mat = opencl_manager.GetDeviceBuffer<GGEMSMaterialTables>(material_tables, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, sizeof(GGEMSMaterialTables), thread_index);

  // Loop over every activated material, build its stopping-power curve from its real elemental composition
  for (GGsize material_index = 0; material_index < cs->number_of_materials_; ++material_index) {
    // Gather this material's elements (Z, mass fraction) from GGEMSMaterialTables
    std::vector<ElementComposition> composition;
    GGsize n_elements = mat->number_of_chemical_elements_[material_index];
    GGsize elem_offset = mat->index_of_chemical_elements_[material_index];
    for (GGsize e = 0; e < n_elements; ++e) {
      ElementComposition ec;
      ec.Z = mat->atomic_number_Z_[elem_offset + e];
      ec.weight_fraction = mat->mass_fraction_[elem_offset + e];
      composition.push_back(ec);
    }

    for (GGsize i = 0; i < cs->number_of_bins_; ++i) {
      GGfloat energy_MeV = cs->energy_bins_[i]; // already stored in GGEMS internal energy units (MeV-based)
      GGfloat sp = ComputeBetheStoppingPower(energy_MeV, composition);
      cs->electron_stopping_power_[material_index * MAX_CROSS_SECTION_TABLE_NUMBER_BINS + i] = sp;
    }
  }

  opencl_manager.ReleaseDeviceBuffer(material_tables, mat, thread_index);
  opencl_manager.ReleaseDeviceBuffer(particle_cross_sections, cs, thread_index);

  GGcout("GGEMSElectronStoppingPower", "BuildCrossSectionTables", 1) << "Electron stopping power tables built for " << cs->number_of_materials_ << " material(s)" << GGendl;
}

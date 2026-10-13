#ifndef GUARD_GGEMS_NAVIGATORS_GGEMSELECTRONNAVIGATOR_HH
#define GUARD_GGEMS_NAVIGATORS_GGEMSELECTRONNAVIGATOR_HH

// ************************************************************************
// * This file is part of GGEMS.                                          *
// * Extension: Electron navigator (CSDA straight-line stepping, v1)      *
// * Added as part of thesis work extending GGEMS to electron transport   *
// ************************************************************************

/*!
  \file GGEMSElectronNavigator.hh

  \brief Functions for electron navigation (CSDA stepping), only for OpenCL kernel usage.
  v1 simplification: straight-line stepping, no multiple scattering, no bremsstrahlung.
*/

#ifdef __OPENCL_C_VERSION__

#include "GGEMS/physics/GGEMSPrimaryParticles.hh"
#include "GGEMS/physics/GGEMSParticleCrossSections.hh"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*!
  \fn inline GGfloat GetElectronStoppingPower(global GGEMSParticleCrossSections const* particle_cross_sections, GGuchar const index_material, GGint const particle_id, global GGEMSPrimaryParticles* primary_particle)
  \brief Look up the electron total stopping power (MeV.cm2/g) at the particle's current energy, for a given material,
  by binary-searching the same log-scale energy_bins_ used for photon cross-sections.
*/
inline GGfloat GetElectronStoppingPower(
  global GGEMSParticleCrossSections const* particle_cross_sections,
  GGuchar const index_material,
  GGint const particle_id,
  global GGEMSPrimaryParticles* primary_particle)
{
  GGint energy_id = BinarySearch(primary_particle->E_[particle_id], particle_cross_sections->energy_bins_, particle_cross_sections->number_of_bins_);
  return particle_cross_sections->electron_stopping_power_[index_material * MAX_CROSS_SECTION_TABLE_NUMBER_BINS + energy_id];
}

#endif

#endif // GUARD_GGEMS_NAVIGATORS_GGEMSELECTRONNAVIGATOR_HH

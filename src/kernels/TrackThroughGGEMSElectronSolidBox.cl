// ************************************************************************
// * This file is part of GGEMS.                                          *
// * Extension: Electron transport kernel (CSDA straight-line stepping)   *
// * Added as part of thesis work extending GGEMS to electron transport   *
// * v1 LIMITATIONS: no multiple scattering, no bremsstrahlung, straight  *
// * line stepping only, fixed number of steps across CSDA range          *
// ************************************************************************

/*!
  \file TrackThroughGGEMSElectronSolidBox.cl

  \brief OpenCL kernel tracking electrons within solid box using CSDA (continuous slowing down approximation)
*/

#include "GGEMS/physics/GGEMSPrimaryParticles.hh"
#include "GGEMS/geometries/GGEMSSolidBoxData.hh"
#include "GGEMS/geometries/GGEMSRayTracing.hh"
#include "GGEMS/materials/GGEMSMaterialTables.hh"
#include "GGEMS/physics/GGEMSParticleCrossSections.hh"
#include "GGEMS/randoms/GGEMSRandom.hh"
#include "GGEMS/maths/GGEMSMatrixOperations.hh"
#include "GGEMS/navigators/GGEMSElectronNavigator.hh"

/*!
  \fn kernel void track_through_ggems_electron_solid_box(...)
  \brief OpenCL kernel tracking electrons within solid box, CSDA straight-line approximation (v1)
*/
kernel void track_through_ggems_electron_solid_box(
  GGsize const particle_id_limit,
  global GGEMSPrimaryParticles* primary_particle,
  global GGEMSSolidBoxData const* solid_box_data,
  global GGuchar const* label_data,
  global GGEMSParticleCrossSections const* particle_cross_sections,
  global GGEMSMaterialTables const* materials,
  GGfloat const threshold,
  GGfloat const step_length_mm
  #ifdef HISTOGRAM
  ,global GGint* histogram
  #endif
)
{
  GGsize global_id = get_global_id(0);
  if (global_id >= particle_id_limit) return;
  if (primary_particle->solid_id_[global_id] != solid_box_data->solid_id_) return;
  if (primary_particle->pname_[global_id] != ELECTRON) return; // only electrons handled by this kernel
  if (primary_particle->status_[global_id] == DEAD) return;

  GGfloat3 global_position = {primary_particle->px_[global_id], primary_particle->py_[global_id], primary_particle->pz_[global_id]};
  GGfloat3 global_direction = {primary_particle->dx_[global_id], primary_particle->dy_[global_id], primary_particle->dz_[global_id]};
  GGfloat3 local_position = GlobalToLocalPosition(&solid_box_data->obb_geometry_.matrix_transformation_, &global_position);
  GGfloat3 local_direction = GlobalToLocalDirection(&solid_box_data->obb_geometry_.matrix_transformation_, &global_direction);

  GGfloat3 border_min = solid_box_data->obb_geometry_.border_min_xyz_;
  GGfloat3 border_max = solid_box_data->obb_geometry_.border_max_xyz_;

  // Material density (v1: material index 0, single-material box, matches current stopping-power table limitation)
  GGuchar material_id = 0;
  GGfloat density_g_cm3 = materials->density_of_material_[material_id]; // already stored in g/cm3, no conversion needed

  // CSDA straight-line stepping loop
  do {
    // Energy-loss this step: stopping_power [MeV.cm2/g] * density [g/cm3] * step_length [cm]
    GGfloat stopping_power = GetElectronStoppingPower(particle_cross_sections, material_id, global_id, primary_particle);
    GGfloat step_length_cm = step_length_mm / 10.0f;
    GGfloat energy_loss = stopping_power * density_g_cm3 * step_length_cm * MeV;

    // Don't lose more energy than the particle has
    GGfloat current_energy = primary_particle->E_[global_id];
    if (energy_loss > current_energy) energy_loss = current_energy;

    // Move particle in straight line (v1: no multiple-scattering deflection)
    local_position = local_position + local_direction * step_length_mm;

    // Deposit the lost energy locally (v1 simplification: no delta-ray/secondary transport,
    // energy is deposited at the midpoint of this step -- a coarse local-deposition model,
    // same limitation photon Compton/Photoelectric currently have)
    #ifdef HISTOGRAM
    if (histogram) {
      // (voxel index calculation would go here, matching the photon kernel's histogram pattern)
    }
    #endif

    primary_particle->E_[global_id] = current_energy - energy_loss;

    // Check if out of box
    if (!IsParticleInAABB(&local_position, border_min.x, border_max.x, border_min.y, border_max.y, border_min.z, border_max.z, GEOMETRY_TOLERANCE)) {
      primary_particle->solid_id_[global_id] = -1;
      break;
    }

    // Check energy threshold -- deposit remaining energy and kill particle
    if (primary_particle->E_[global_id] <= threshold) {
      primary_particle->E_[global_id] = 0.0f;
      primary_particle->status_[global_id] = DEAD;
    }

  } while (primary_particle->status_[global_id] == ALIVE);

  // Store final position/direction back in global frame
  global_position = LocalToGlobalPosition(&solid_box_data->obb_geometry_.matrix_transformation_, &local_position);
  primary_particle->px_[global_id] = global_position.x;
  primary_particle->py_[global_id] = global_position.y;
  primary_particle->pz_[global_id] = global_position.z;
}

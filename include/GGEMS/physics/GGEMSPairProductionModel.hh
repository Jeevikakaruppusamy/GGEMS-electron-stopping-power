#ifndef GUARD_GGEMS_PHYSICS_GGEMSPAIRPRODUCTIONMODEL_HH
#define GUARD_GGEMS_PHYSICS_GGEMSPAIRPRODUCTIONMODEL_HH

#ifdef __OPENCL_C_VERSION__

inline void PairProductionSampleSecondaries(
  global GGEMSPrimaryParticles* primary_particle,
  GGint const particle_id
)
{
  primary_particle->status_[particle_id] = DEAD;
  primary_particle->E_[particle_id] = 0.0f;
}

#endif

#endif // GUARD_GGEMS_PHYSICS_GGEMSPAIRPRODUCTIONMODEL_HH

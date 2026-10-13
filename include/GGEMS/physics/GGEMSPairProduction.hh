#ifndef GUARD_GGEMS_PHYSICS_GGEMSPAIRPRODUCTION_HH
#define GUARD_GGEMS_PHYSICS_GGEMSPAIRPRODUCTION_HH

#include "GGEMS/physics/GGEMSEMProcess.hh"

class GGEMS_EXPORT GGEMSPairProduction : public GGEMSEMProcess
{
  public:
    GGEMSPairProduction(std::string const& primary_particle, bool const& is_secondary);
    ~GGEMSPairProduction(void);

    GGEMSPairProduction(GGEMSPairProduction const& pair_production) = delete;
    GGEMSPairProduction& operator=(GGEMSPairProduction const& pair_production) = delete;
    GGEMSPairProduction(GGEMSPairProduction const&& pair_production) = delete;
    GGEMSPairProduction& operator=(GGEMSPairProduction const&& pair_production) = delete;

    GGfloat ComputeCrossSectionPerAtom(GGfloat const& energy, GGuchar const& atomic_number) const override;
};

#endif // GUARD_GGEMS_PHYSICS_GGEMSPAIRPRODUCTION_HH

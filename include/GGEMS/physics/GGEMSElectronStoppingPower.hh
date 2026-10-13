#ifndef GUARD_GGEMS_PHYSICS_GGEMSELECTRONSTOPPINGPOWER_HH
#define GUARD_GGEMS_PHYSICS_GGEMSELECTRONSTOPPINGPOWER_HH

// ************************************************************************
// * This file is part of GGEMS.                                          *
// * Extension: Electron stopping-power table (NIST ESTAR based)          *
// * Added as part of thesis work extending GGEMS to electron transport   *
// ************************************************************************

/*!
  \file GGEMSElectronStoppingPower.hh

  \brief GGEMS class handling electron stopping-power tables (from NIST ESTAR data)
*/

/// \cond
#include <vector>
#include <string>
/// \endcond

#include "GGEMS/physics/GGEMSEMProcess.hh"

/*!
  \class GGEMSElectronStoppingPower
  \brief GGEMS class handling electron stopping-power tables, loaded from NIST ESTAR reference data
  instead of computed from an analytic per-atom cross-section formula.
*/
class GGEMS_EXPORT GGEMSElectronStoppingPower : public GGEMSEMProcess
{
  public:
    /*!
      \brief GGEMSElectronStoppingPower constructor
    */
    GGEMSElectronStoppingPower(void);

    /*!
      \brief GGEMSElectronStoppingPower destructor
    */
    ~GGEMSElectronStoppingPower(void);

    GGEMSElectronStoppingPower(GGEMSElectronStoppingPower const& sp) = delete;
    GGEMSElectronStoppingPower& operator=(GGEMSElectronStoppingPower const& sp) = delete;
    GGEMSElectronStoppingPower(GGEMSElectronStoppingPower const&& sp) = delete;
    GGEMSElectronStoppingPower& operator=(GGEMSElectronStoppingPower const&& sp) = delete;

    /*!
      \fn void LoadStoppingPowerTable(std::string const& material_name, std::string const& file_path)
      \param material_name - name of the material this table applies to (e.g. "Water")
      \param file_path - path to the NIST ESTAR-format data file
      \brief Load a NIST ESTAR stopping-power table from file into host memory
    */
    void LoadStoppingPowerTable(std::string const& material_name, std::string const& file_path);

    /*!
      \fn void BuildCrossSectionTables(cl::Buffer* particle_cross_sections, cl::Buffer* material_tables, GGsize const& thread_index)
      \brief Override: fill electron_stopping_power_ table using log-log interpolation of loaded NIST data,
      instead of the per-atom analytic formula used by photon processes.
    */
    void BuildCrossSectionTables(cl::Buffer* particle_cross_sections, cl::Buffer* material_tables, GGsize const& thread_index) override;

  protected:
    /*!
      \brief Not used for table-based stopping power (required override from GGEMSEMProcess base class)
    */
    GGfloat ComputeCrossSectionPerAtom(GGfloat const& energy, GGuchar const& atomic_number) const override;

  private:
    /*!
      \fn GGfloat InterpolateStoppingPower(GGfloat const& energy_MeV) const
      \param energy_MeV - energy in MeV to interpolate at
      \return interpolated total stopping power in MeV.cm2/g
      \brief Log-log linear interpolation of the loaded NIST ESTAR table
    */
    GGfloat InterpolateStoppingPower(GGfloat const& energy_MeV) const;

  public:
    /*!
      \struct ElementComposition
      \brief One element's weight fraction within a compound/mixture material
    */
    struct ElementComposition {
      GGuchar Z;   /*!< Atomic number */
      GGfloat weight_fraction; /*!< Mass fraction of this element in the material */
    };

    /*!
      \fn GGfloat ComputeBetheStoppingPower(GGfloat const& energy_MeV, std::vector<ElementComposition> const& composition) const
      \param energy_MeV - electron kinetic energy in MeV
      \param composition - list of elements and their mass fractions making up the material
      \return total (collision) mass stopping power in MeV.cm2/g, via the Bethe formula and Bragg additivity
      \brief Compute stopping power directly from elemental composition (Bethe formula, ICRU-37 style),
      instead of requiring a per-material NIST ESTAR table. This generalizes to any material whose
      elemental composition is known (e.g. from GGEMS materials.txt), such as the 53 ICRP-110 materials.
      v1 limitation: collision stopping power only (no radiative/bremsstrahlung term), no density effect
      correction -- reasonable below a few MeV, which covers most diagnostic/therapeutic photon-secondary
      electron energies.
    */
    GGfloat ComputeBetheStoppingPower(GGfloat const& energy_MeV, std::vector<ElementComposition> const& composition) const;

  private:
    /*!
      \fn GGfloat BetheStoppingPowerPerElement(GGfloat const& energy_MeV, GGint const& Z, GGfloat const& A, GGfloat const& I_eV) const
      \brief Bethe collision stopping power (MeV.cm2/g) for a single element at a given electron energy
    */
    GGfloat BetheStoppingPowerPerElement(GGfloat const& energy_MeV, GGint const& Z, GGfloat const& A, GGfloat const& I_eV) const;

  private:
    std::vector<GGfloat> nist_energy_MeV_; /*!< Energies from the loaded NIST table, in MeV */
    std::vector<GGfloat> nist_stopping_power_; /*!< Total stopping power (Stot/rho) values, MeV.cm2/g */
    std::string loaded_material_name_; /*!< Name of material the loaded table applies to */
};

#endif // GUARD_GGEMS_PHYSICS_GGEMSELECTRONSTOPPINGPOWER_HH

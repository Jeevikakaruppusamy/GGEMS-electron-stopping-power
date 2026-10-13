import ctypes
from ggems import *

GGEMSVerbosity(0)
opencl_manager = GGEMSOpenCLManager()
materials_database_manager = GGEMSMaterialsDatabaseManager()
processes_manager = GGEMSProcessesManager()
opencl_manager.set_device_to_activate('cpu')
materials_database_manager.set_materials('data/materials.txt')

materials = GGEMSMaterials()
materials.add_material('Water')
materials.initialize()

processes_manager.set_cross_section_table_number_of_bins(220)
processes_manager.set_cross_section_table_energy_min(1.0, 'keV')
processes_manager.set_cross_section_table_energy_max(10.0, 'MeV')

cross_sections = GGEMSCrossSections(materials)
cross_sections.add_process('ElectronIonisation', 'e-')
cross_sections.initialize()

# CSDA-style stepping: start at 1.0 MeV, step down in small energy decrements,
# accumulate path length = sum of (delta_E / stopping_power) -- density in g/cm3 = 1.0 for water
start_energy = 1.0  # MeV
density = 1.0  # g/cm3, water
n_steps = 2000
dE = start_energy / n_steps

energy = start_energy
path_length_cm = 0.0
for i in range(n_steps):
    if energy <= 0.01:  # stop near table's practical lower bound
        break
    sp = cross_sections.get_electron_stopping_power('Water', energy, 'MeV')  # MeV.cm2/g
    if sp <= 0:
        break
    step_cm = dE / (sp * density)
    path_length_cm += step_cm
    energy -= dE

print(f"Computed CSDA-style range for 1.0 MeV electron in water: {path_length_cm:.4f} cm")
print(f"NIST reference CSDA range for 1.0 MeV in water: 0.4384 cm (g/cm2, density=1 so same in cm)")
print(f"Difference: {abs(path_length_cm - 0.4384)/0.4384*100:.1f}%")

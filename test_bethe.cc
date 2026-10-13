#include <iostream>
#include <iomanip>
#include <vector>
#include "GGEMS/physics/GGEMSElectronStoppingPower.hh"

int main() {
    GGEMSElectronStoppingPower sp;

    // Water composition (from GGEMS materials.txt: H f=0.111, O f=0.889)
    std::vector<GGEMSElectronStoppingPower::ElementComposition> water = {
        {"H", 0.111f},
        {"O", 0.889f}
    };

    // Bone composition (from materials_icrp110.txt)
    std::vector<GGEMSElectronStoppingPower::ElementComposition> bone = {
        {"H", 0.036f}, {"C", 0.159f}, {"N", 0.042f}, {"O", 0.448f},
        {"Na", 0.003f}, {"Mg", 0.002f}, {"P", 0.094f}, {"S", 0.003f}, {"Ca", 0.213f}
    };

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "=== WATER: Bethe (collision only) vs NIST (total, incl. radiative) ===\n";
    std::cout << "Energy(MeV) | Bethe (collision) | NIST total (reference)\n";
    float energies[] = {0.1f, 0.5f, 1.0f, 2.0f, 10.0f};
    float nist_total[] = {4.097f, 2.032f, 1.858f, 1.848f, 2.148f};
    for (int i = 0; i < 5; i++) {
        float bethe_val = sp.ComputeBetheStoppingPower(energies[i], water);
        std::cout << energies[i] << "        | " << bethe_val << "             | " << nist_total[i] << "\n";
    }

    std::cout << "\n=== BONE (ICRP-110 composition): Bethe collision stopping power ===\n";
    std::cout << "Energy(MeV) | Stopping power (MeV.cm2/g)\n";
    for (int i = 0; i < 5; i++) {
        float bethe_val = sp.ComputeBetheStoppingPower(energies[i], bone);
        std::cout << energies[i] << "        | " << bethe_val << "\n";
    }

    return 0;
}

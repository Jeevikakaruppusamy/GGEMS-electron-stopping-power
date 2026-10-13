#include <iostream>
#include <iomanip>
#include "GGEMS/physics/GGEMSElectronStoppingPower.hh"

int main() {
    GGEMSElectronStoppingPower sp;
    sp.LoadStoppingPowerTable("Water", "/home/cwr/GGEMS/data/electron_stopping_power/water_estar.txt");

    // Known exact values from the NIST ESTAR table (Table A.3, water):
    struct TestPoint { float energy; float expected; };
    TestPoint points[] = {
        {0.1f, 4.097f},
        {0.5f, 2.032f},
        {1.0f, 1.858f},
        {2.0f, 1.848f},
        {10.0f, 2.148f}
    };

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Energy(MeV) | Interpolated | Expected | Match?\n";
    std::cout << "------------------------------------------------\n";
    for (auto& p : points) {
        float result = sp.InterpolateStoppingPower(p.energy);
        bool match = std::abs(result - p.expected) < 0.01f;
        std::cout << p.energy << "       | " << result << "       | " << p.expected << "   | " << (match ? "YES" : "NO") << "\n";
    }

    // Also test a value BETWEEN table points, to verify interpolation (not just exact lookup)
    float mid_energy = 0.3f; // between 0.2 (2.784) and 0.4 (2.144) in the table -> not exact, real interpolation test
    std::cout << "\nInterpolation-only test at " << mid_energy << " MeV: " << sp.InterpolateStoppingPower(mid_energy) << " MeV.cm2/g\n";
    std::cout << "(Should land between neighboring table values, and match table's own 0.3 MeV entry: 2.349)\n";

    return 0;
}

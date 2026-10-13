#include <iostream>
#include <fstream>
#include <cmath>
#include "GGEMS/physics/GGEMSElectronStoppingPower.hh"

int main() {
    GGEMSElectronStoppingPower sp;
    sp.LoadStoppingPowerTable("Water", "/home/cwr/GGEMS/data/electron_stopping_power/water_estar.txt");

    std::ofstream out("stopping_power_output.csv");
    out << "energy_MeV,stopping_power\n";

    // Sample densely in log space from 0.001 to 10 MeV, using GGEMS's own
    // interpolation logic (the same InterpolateStoppingPower the kernel will use)
    float e_min = 0.001f, e_max = 10.0f;
    int n_points = 200;
    for (int i = 0; i < n_points; i++) {
        float log_e = logf(e_min) + (logf(e_max) - logf(e_min)) * i / (n_points - 1);
        float energy = expf(log_e);
        float value = sp.InterpolateStoppingPower(energy);
        out << energy << "," << value << "\n";
    }
    out.close();

    std::cout << "Wrote " << n_points << " points to stopping_power_output.csv" << std::endl;
    return 0;
}

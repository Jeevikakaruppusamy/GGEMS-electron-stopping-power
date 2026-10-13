# GGEMS Electron Stopping Power

Experimental tests and examples for electron stopping-power and CSDA-range work with [GGEMS](https://github.com/GGEMS/GGEMS).

## Files in this repository

- `test_sp.cc` — checks stopping-power table lookups against selected NIST ESTAR water values and exercises interpolation.
- `test_bethe.cc` — compares Bethe collision stopping-power calculations for water and bone with reference values.
- `test_csda_range.py` — estimates the CSDA-style range for a 1 MeV electron in water using GGEMS cross sections.
- `plot_sp.cc` — writes interpolated stopping-power values to `stopping_power_output.csv` for plotting.

## Requirements and setup

These are companion tests, not a standalone GGEMS build. Running the C++ examples requires a GGEMS source/build environment with its headers, implementation, and dependencies available. The Python example requires the GGEMS Python module, OpenCL support, and the GGEMS materials data file.

The C++ examples currently refer to the stopping-power data file at `/home/cwr/GGEMS/data/electron_stopping_power/water_estar.txt`. Update that path to the location of your GGEMS checkout before running them. The Python example expects `data/materials.txt` to be available from its working directory.

## Notes

The values in these examples are exploratory comparisons, not a guarantee of validated physics results. Rebuild executables from their source in your own GGEMS environment rather than relying on uploaded binaries.

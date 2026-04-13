# WellRTSim

WellRTSim is a MATLAB wellbore reactive transport simulator for coupled flow, heat transfer, and chemistry calculations.

## Quick Start

1. Open `main.m` in MATLAB.
2. Set `SimDir` in `main.m` to the simulation case you want to run.
3. Verify the PHREEQC database path in the selected case `chemistry.md` file.
4. Run `main.m`.

Simulation inputs live under `Simulation/`, results are written to the selected case directory, and the plotting/viewer utilities are in `viewer/`.

## Notes

- The default configured case in `main.m` is `Simulation/Krafla/KJ-9/`.
- `user_manual.pdf` contains additional project documentation.
- `PhreeqcMatlab/` provides the MATLAB interface used for chemistry calculations.

## License

WellRTSim is licensed under CC BY-NC-SA 4.0 for academic and research use.

For commercial licensing inquiries, contact `oleg.melnik@earth.ox.ac.uk` or `oemelnik@gmail.com`.

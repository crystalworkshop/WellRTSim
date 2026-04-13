# WellRTSim

WellRTSim is a MATLAB wellbore reactive transport simulator for coupled flow, heat transfer, and chemistry calculations.

## Quick Start

1. Install the external `PhreeqcMatlab` package from <https://github.com/simulkade/PhreeqcMatlab>.
2. In MATLAB, run the `startup.m` file from that `PhreeqcMatlab` checkout so `IPhreeqc` and the PHREEQC databases are added to the MATLAB path.
3. Open `main.m` in MATLAB.
4. Set `SimDir` in `main.m` to the simulation case you want to run.
5. Verify the PHREEQC database entry in the selected case `chemistry.md` file.
6. Run `main.m`.

Simulation inputs live under `Simulation/`, results are written to the selected case directory, and the plotting/viewer utilities are in `viewer/`.

## Notes

- The default configured case in `main.m` is `Simulation/Krafla/KJ-9/`.
- The default chemistry setup uses `phreeqc.dat`, which WellRTSim resolves from the active MATLAB path after `PhreeqcMatlab/startup.m` has been run.
- `user_manual.pdf` contains additional project documentation.

## License

WellRTSim is licensed under CC BY-NC-SA 4.0 for academic and research use.

For commercial licensing inquiries, contact `oleg.melnik@earth.ox.ac.uk` or `oemelnik@gmail.com`.

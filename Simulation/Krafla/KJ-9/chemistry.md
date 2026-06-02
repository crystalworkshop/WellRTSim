# Chemistry Setup (Krafla KJ-9)

## Aqueous species

                   Ca, C(4), Si,  Na,   K,   S(6), S(-2),   Cl
concentration_ppm: 22, 400, 409, 139, 17.5, 305,  64.1   25.6
Suffix

    "", "as HCO3+", "", "", "", "as SO4", "as H2S", ""

## Gases:

CO2(g), H2S(g)
partition: 2 # 1 - phreeqc, 2 - defined in matlab

## Minerals

Kinetics: Calcite, Quartz, Anhydrite
formula: CaCO3, SiO2, CaSO4
tau_s: 8e4, 1e6, 1e6
density_kg_m3: 2710, 2650, 2980

## Template file

chemistry.pht
generate phreeqc: true

## Database file path
/usr/local/share/doc/iphreeqc/database/phreeqc.dat
# llnl.dat
# /absolute/path/to/PhreeqcMatlab/database/phreeqc.dat

# Simulation Parameters
# Case: Tonkin "Problem 3 / East Mesa 6-1" - 198.5 degC low-enthalpy production, WITH heat loss to formation.
# Same IC/parameters as LO_T (LT_no_heat_loss) except reservoir_heat_flux = true.
# Extracted from LT_heat_loss.json (boundaries, feedzone, time) and LT_mesh.json (geometry).

## Grid Parameters

n = 64 # Number of grid points (matches mesh NEL = 50)
maxiter = 20 # Maximum number of iterations
epsQ = 1e-6 # Convergence tolerance
P_unit = bar # bar or MPa

## Well Steady-state Conditions

P_top = 2.057180 # Outlet pressure, bar (final top P from LT_no_heat_loss.h5 @ 1600 d)
iBC_top = 3 # 1 - WHP.csv, 2 - specified flow rate, 3 - fixed outlet pressure. JSON wellhead = production/mass_flux -> 2.
whp_flow_control = 1 # 1 - relieve WHP toward atmospheric when top flow goes negative (floored at P_atm)
whp_ctrl_relax = 0.1 # fraction of (WHP - P_atm) shed per timestep while backflow persists
feed = 0  # 0 - PI based feed (JSON feedzone type = productivity), 2 - fixed bottom P/H, 3 - specified feedzone flow.
IC_switch =3  # 1 - top-down steady init, 2 - bottom-up steady init, 3 - prescribed InitPressure/InitTemperature profiles
P_bot = 200 # 95.882210 # Bottom pressure, bar (feedzone reservoir P = 21 MPa, closed bottomhole)
T_bot = 198.5  # Bottom temperature, degC (feedzone fluid temperature)
H_top = 854 # Top enthalpy guess, kJ/kg (liquid at 198.5 degC, 21 MPa)
Q_init = 0.0 # target flow rate kg/s (wellhead mass_flux plateau after soft-start)

## Fluid Chemistry

calc_chem = 0 # 1 - calculate chemical transport (disabled for this case)

## Physical Properties

g = 9.81 # Gravitational acceleration (JSON g = -9.81)
rho_r = 2700 # Rock density
C_r = 1000 # Rock heat capacity
k_r = 2.42 # Rock thermal conductivity W/m/K
H_q = 62.6 # heat transfer coefficient U [W/m^2/K], heat loss ON. U = k_c/(r*ln(r_c/r)) with cement k_c=2.42 W/m/K, r=0.11075 m (inner radius), r_c=0.157 m (casing OD/2). U = 25.88*k_c.

## Time Settings

tunit = h     # time unit d - days, h - hours, s - seconds
tfin = 24 # Final time in time units (JSON time.stop = 138240000 s)
dt_max = 180  # Max time step in seconds (JSON time.step.maximum.size)
dt = 1.0 # initial time step, s (JSON time.step.initial)
t_adjust = 100 # no heat transfer during adjustment, s
dt_increment = 1.15 # increment in time step
pltf = 10    # Plot frequency
t_save = 0.0525   # save frequency in tunit
save_csv = 0 # 1 save csv profiles 0 - only results.h5
results_prefix = LT_

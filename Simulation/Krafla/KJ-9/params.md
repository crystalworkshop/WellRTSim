# Simulation Parameters

## Grid Parameters

n = 128 # Number of grid points
maxiter = 10 # Maximum number of iterations
epsQ = 1e-4 # Convergence tolerance
P_unit = bar # bar or MPa

## Well Steady-state Conditions

P_top = 5.52 # Outlet pressure, bar
iBC_top = 3 # 1 - specified wellhead pressure from WHP.csv, 2 - specified flow rate, 3 - fixed outlet pressure equal to P_top.
feed = 0  # 0 - PI based feed, 2 - fixed bottom pressure/enthalpy, 3 - specified feedzone flow rate/enthalpy.
IC_switch =2  # 1 - top-down steady init, 2 - bottom-up steady init
P_bot = 77.6 # Bottom pressure, P_unit (used in bottom-up initialization)
T_bot = 211.  # Bottom temperature, degC (used in bottom-up initialization)
H_top = 1510 # Top enthalpy, kJ/kg
Q_init = 8. # target flow rate kg/s.

## Fluid Chemistry

calc_chem = 1 # 1- calculate chemical transport
stat_chem  = 600  # time in sec to start chemical transport to allow IC equilibration
chem_source = 1 # 0 - disable lateral chemical source from feedzones, 1 - enable
chem_semi_implicit_max_iter = 3 # internal chemistry transport-PHREEQC iterations per outer timestep
chem_semi_implicit_relaxation = 0.35 # relaxation factor for partition-coefficient update, 0..1
eps_scale = 2e-3 # scale roughness, m
chemistry setup = chemistry.md

## Physical Properties

g = 9.81 # Gravitational acceleration
rho_r = 2700 # Rock density
C_r = 1000 # Rock heat capacity
k_r = 2.42 # Rock thermal conductivity W/m/K
H_q =3e2 # heat conduction coefficient

## Time Settings

tunit = d     #time unit d - days, h - hours, s - seconds
tfin = 180 # Final time in time unints
dt_max =  36000   # Max time step in seconds
dt = 0.5 #initial time step, s
t_adjust =100 # no heat transfer during ajustment, s
dt_increment = 1.15 #increment in time step
pltf = 10    # Plot frequency
t_save = 5   # save frequency in tunit
save_csv = 0 # 1 save csv profiles 0 - only results.h5
results_prefix = Krafla_

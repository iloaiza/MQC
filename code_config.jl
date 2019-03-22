#= This file holds the configurations for the code. Anything can be changed by using
the codewords as seen in constant definitions.

=#
## CONTROL VARIABLES FOR DYNAMICS
const rk_tol = 1e-5 #default tolerance for each step of runge-kutta 45 adaptative algorithm
const dt_min = 1e-3 #default minimum timestep for rk45 (2.42e-5fs)
const dt_max = 20.0 #default maximum timestep for rk45 (0.12fs)

## VARIABLES FOR DEBUGGING AND PRINTING CHECKS THROUGHOUT DYNAMICS
const print_phase = false #make true to print warnings when phase (from dia->adi transformation) is being used as sign and not just rounding
const tol = 1e-3 #default tolerance for sanity checks (energy and norm conservation). 1e-3 default -> less than 0.1% variation. Activate with sanity_checks
const sanity_checks = true #will perform sanity checks by the end of every trajectory, comparing energy change and norm change to the defined tol
#NOT IMPLEMENTED# const high_verbose_sanity_checks = false #turn true for sanity check every timestep or checkpoint (for tracking region where dynamics break)
const time_print = true #shows the maximum, minimum, and mean timestep per trajectory in output

## FOR AUTOMATIC OUTPUT PLOTTING
const plot_out = false #turn true for automatic generation of default plots

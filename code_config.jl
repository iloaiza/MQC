#= This file holds the configurations for the code. Anything can be changed by using
the codewords as seen in constant definitions after the include section of any Initial_data file.

=#
## CONTROL VARIABLES FOR DYNAMICS
const rk_tol = 1e-5 #default tolerance for each step of runge-kutta 45 adaptative algorithm
const dt_min = 1e-3 #default minimum timestep for rk45 (2.42e-5fs)
const dt_max = 10.0 #default maximum timestep for rk45 (0.12fs)

## VARIABLES FOR DEBUGGING AND PRINTING CHECKS THROUGHOUT DYNAMICS
const print_phase = false #make true to print warnings when phase (from dia->adi transformation) is being used as sign and not just rounding
const tol = 1e-3 #default tolerance for sanity checks (energy and norm conservation). 1e-3 default -> less than 0.1% variation. Activate with sanity_checks
const sanity_checks = true #will perform sanity checks by the end of every trajectory, comparing energy change and norm change to the defined tol
#NOT IMPLEMENTED# const high_verbose_sanity_checks = false #turn true for sanity check every timestep or checkpoint (for tracking region where dynamics break)
const time_print = true #shows the maximum, minimum, and mean timestep per trajectory in output

## WALLTIMES AND RESOLUTION OF SAVES
const walltime = 10000 #walltime, in a.u. (~242fs) for simulations which continue until the trajectory has left interaction zone
const checkpoints = 100 #default number of checkpoints for DYNAMICS simulations
FIRST_RUN = false #true to create savefile with initial energy (also creates folder for K_SIMULATIONS)

## FOR AUTOMATIC OUTPUT PLOTTING
const plot_out = true #turn true for automatic generation of default plots
const plot_method = "plotlyjs" #for plotting method (in Plots package)
const HISTO_RES = 150 #resolution of histograms (i.e. number of bars)
const SH_eval = true #when doing SH statistics, chooses final adiabatic state (false for mean of electronic properties)

#= This file holds the configurations for the code. Anything can be changed by using
the codewords as seen in constant definitions after the include section of any Initial_data file.

=#
## DYNAMICAL REGIMES
initial_dist = false; #default initial dist, false corresponds to constant_dist (all initial positions and momenta as R0 and p0)

## CONTROL VARIABLES FOR ADAPTATIVE DYNAMICS (RK54 (way better results than 45!))
const RelTol = [1e-5,1e-5,1e-5,1e-5] #default relative tolerance for each step of runge-kutta 45 adaptative algorithm
##### RelTol = [R_tol,p_tol,C_tol,mem_tol], where R is nuc_position, p is momentum, C is electronic coeffs and mem is memory for friction integrations
const AbsTol = [1e-6,1e-6,1e-6,1e-6] #default absolute tolerances for rk54 adaptative timestep algorithm, follows RelTol conventions
const dt_min = 1e-4 #default minimum timestep for rk54 (1e-4 a.u. = 2.42e-6 fs)
const dt_max = 1.0 #default maximum timestep for rk54 (1 a.u. = 0.024 fs)

## VARIABLES FOR DEBUGGING AND PRINTING CHECKS THROUGHOUT DYNAMICS
const print_phase = true #make true to print warnings when phase (from dia->adi transformation) is being used as sign and not just rounding
const tol = 1e-3 #default tolerance for sanity checks (energy and norm conservation). 1e-3 default -> less than 0.1% variation. Activate with sanity_checks
const sanity_checks = true #will perform sanity checks by the end of every trajectory, comparing energy change and norm change to the defined tol
const high_verbose = false #turn true for sanity check every timestep or checkpoint (for tracking region where dynamics break). Use in conjunction with time_print
const sanity_breaks = false #true: stops when sanity check is broken, false to just display warning

## RUNGE-KUTTA OPTIONS
const time_print = false #true: shows the maximum, minimum, and mean timestep per trajectory in output (for adaptative timestep!)
const fixed_step = false #true for fixing timestep, removing adaptative timestep and using dt from initial file (rk5)
const adaptative_verbose = false #true: prints calculated errors and trial step every runge45 run (only for adaptative timestep!)
const rk5 = false #when using fixed step, turn true to use rk5 using the coefficients from rk54 (uses rk4 for error control and integrates with the more exact rk5 algorithm)
const rk4 = true #when using fixed step, turn true to use rk4


## WALLTIMES AND RESOLUTION OF SAVES
const walltime = 20000 #walltime, in a.u. (~484fs) for simulations which continue until the trajectory has left interaction zone
const checkpoints = 100 #default number of checkpoints for DYNAMICS simulations
FIRST_RUN = true #true to create savefile with initial energy (also creates folder for K_SIMULATIONS), false to save on top of already created files

## FOR AUTOMATIC OUTPUT PLOTTING
const plot_out = true #turn true for automatic generation of default plots
const plot_method = "gr" #for plotting method (in Plots package)
const HISTO_RES = 150 #resolution of histograms (i.e. number of bars)
const SH_eval = true #when doing SH statistics, chooses final adiabatic state (false for mean of electronic properties)

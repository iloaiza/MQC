#= This file holds the configurations for the code. Anything can be changed by using
the codewords as seen in constant definitions after the include section of any Initial_data file.

=#
## DYNAMICAL REGIMES
initial_dist = false; #default initial dist, false corresponds to constant_dist (all initial positions and momenta as R0 and p0)

## CONTROL VARIABLES FOR ADAPTATIVE DYNAMICS (RK54 (way better results than 45!))
const RelTol = [1e-5,1e-5,1e-5,1e-5] #default relative tolerance for each step of adaptative algorithm
##### RelTol = [R_tol,p_tol,C_tol,mem_tol], where R is nuc_position, p is momentum, C is electronic coeffs and mem is memory for friction integrations
const AbsTol = [1e-6,1e-6,1e-6,1e-6] #default absolute tolerances for adaptative timestep algorithm,s follows RelTol conventions
const dt_min = 1e-4 #default minimum timestep for rk54 (1e-4 a.u. = 2.42e-6 fs)
const dt_max = 1.0 #default maximum timestep for rk54 (1 a.u. = 0.024 fs)

## VARIABLES FOR DEBUGGING AND PRINTING CHECKS THROUGHOUT DYNAMICS
const print_phase = true #make true to print warnings when phase (from dia->adi transformation) is being used as sign and not just rounding
const tol = 1e-3 #default tolerance for sanity checks (energy and norm conservation). 1e-3 default -> less than 0.1% variation. Activate with sanity_checks
const sanity_checks = true #will perform sanity checks by the end of every trajectory, comparing energy change and norm change to the defined tol
const high_verbose = false #turn true for sanity check every timestep or checkpoint (for tracking region where dynamics break). Use in conjunction with time_print
const sanity_breaks = false #true: stops when sanity check is broken, false to just display warning

## INTEGRATION OPTIONS
const time_print = false #true: shows the maximum, minimum, and mean timestep per trajectory in output (for adaptative timestep)
const fixed_step = false #true for fixing timestep, removing adaptative timestep and using dt from initial file (rk5)
const adaptative_verbose = false #true: prints calculated errors and trial step every rungestep run (for adaptative timestep debugging)
const rk5 = false #when using fixed step rk5, adaptative timestep uses Dorand-Prince 5(4)
const rk4 = true #when using fixed step, turn true to use rk45


## WALLTIMES AND RESOLUTION OF SAVES
const walltime = 20000 #walltime, in a.u. (~484fs) for simulations which continue until the trajectory has left interaction zone
if @isdefined(flags)
	const checkpoints = flags #default number of checkpoints for DYNAMICS simulations
else
	const checkpoints = 100 #default number of flags is 100
end
const FIRST_RUN = true #true to move old files and create savefile with initial energy (also creates folder for K_SIMULATIONS), false to save on top of already created files
const IND_SAVES = true #true to save indiviually each trajectory during integration, deletes individual files when simulation complete. false for just final save

## FOR AUTOMATIC OUTPUT PLOTTING
const plot_out = true #turn true for automatic generation of default plots
const plot_method = "pyplot" #for plotting method (in Plots package)
const HISTO_RES = 150 #resolution of histograms (i.e. number of bars)
const HISTO_MINS = [-28,-28]
const HISTO_MAXS = [23,23]
const SH_eval = true #when doing SH statistics, chooses final adiabatic state (false for mean of electronic properties)
const plot_ini = false #true for plotting initial distribution in histograms for DYNAMICS simulations
const DO_DYN = true #false skips dynamics and passes to plots directly
const BO_COMP = false #true shows mean deviation of trajectories with respect to Born-Oppenheimer dynamics
const EL_PLOTS = false #true plots electronic populations vs time
const EL_ES_PLOTS = false #true for plotting all electronic populations, false for just ground state
const HISTO_PLOTS = true #true for plotting nuclear distribution histograms at each flag


## CMFSH tolerance δf for frictional state on-the-fly detection
const δf_CMFSH = 1.5e-4

## Cleaning for removing dynamics from savefile before starting dynamical routine
CLEAN_DYNS = String[]
const tol = 1e-3 #default tolerance for sanity checks (energy and norm conservation). 1e-3 default -> less than 0.1% variation
const print_phase = false #make true to print warnings when phase is being used as sign and not just rounding

const sanity_checks = true #will perform sanity checks by the end of every trajectory, comparing energy change and norm change
#to the defined tol

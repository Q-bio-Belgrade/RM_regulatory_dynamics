import numpy as np
from scipy.stats import iqr
from multiprocessing import Pool, cpu_count

def gillespie_const(tmax, dt, phiR, phiM, n, R0, M0):
    R = R0
    p = int(round(n))
    M = M0

    t = 0.0
    t_save = dt
    output = []
    output.append([t, R, M, p])

    while t < tmax:
        # Propensities
        p_plus = p
        p_minus = p
        R_plus = p * phiR
        R_minus = R
        M_plus = p * phiM
        M_minus = M

        p_t = p_plus + p_minus + R_plus + R_minus + M_plus + M_minus

        if p_t == 0:
            # No more reactions can occur, system is static.
            # Fill remaining times with current state.
            while t_save <= tmax:
                output.append([t_save, R, M, p])
                t_save += dt
            break

        # Time to next reaction
        tau = np.log(1.0/np.random.rand()) / p_t
        t = t + tau

        # Selecting reaction
        r = np.random.rand()
        cumulative = 0.0
        if r < (cumulative := cumulative + p_plus/p_t):
            p += 1
        elif r < (cumulative := cumulative + p_minus/p_t):
            p -= 1
        elif r < (cumulative := cumulative + R_plus/p_t):
            R += 1
        elif r < (cumulative := cumulative + R_minus/p_t):
            R -= 1
        elif r < (cumulative := cumulative + M_plus/p_t):
            M += 1
        else:
            M -= 1

    # If t goes beyond multiple increments of t_save, record all intermediate steps
    while t >= t_save and t_save <= tmax:
        output.append([t_save, R, M, p])
        t_save += dt

    return np.array(output)

constitutive_simulation_params = np.loadtxt("constitutive_simulation_params.csv", delimiter=",")
phiLambdaR , phiLambdaM , n = constitutive_simulation_params

tmax = 12
R0 = 0
M0 = 0
ts_e = 10
N = 100 #number of points recorded in the simulation
RunsDist = 2000
dt = tmax / N
workers = 4 #set number of workers

max_workers = cpu_count()

if not( 1 <= workers <= max_workers ):
    workers = max_cores


def run_sim(_):
    return gillespie_const(tmax, dt, phiLambdaR, phiLambdaM, n, R0, M0)

if __name__ == '__main__':

    # Generate a seed based on the current time
    #seed_value = int(time.time())
    #np.random.seed(seed_value)

    # Save this seed to a file so you can reuse it later
    #with open("saved_seed_constitutive.txt", "w") as f:
        #f.write(str(seed_value))

    with open("saved_seed_constitutive.txt", "r") as f:
        loaded_seed = int(f.read())

    # Set the seed
    np.random.seed(loaded_seed)

    # Parallel execution of simulations
    with Pool(processes=max_cores) as pool:
        Simulation_constitutive = pool.map(run_sim, range(RunsDist))

    MR_vec_e = []

    for i in range(RunsDist):
        output = Simulation_constitutive[i]
        ts = output[:, 0]
        R_constitutive = output[:, 1]
        M_constitutive = output[:, 2]

        # Find first index where ts > ts_e
        indices = np.where(ts > ts_e)[0]
        index_e = indices[0]
        #if len(indices) > 0:
            #index_e = indices[0]
        if R_constitutive[index_e] > 0 and M_constitutive[index_e] > 0:
            MR_vec_e.append(M_constitutive[index_e] / R_constitutive[index_e])

    MR_vec_e = np.array(MR_vec_e)
    MR_Median = np.median(MR_vec_e)
    M_R_iqr = iqr(MR_vec_e)
    CV_MR_constitutive = M_R_iqr / MR_Median
    print(CV_MR_constitutive)

    # Save the parameters in CSV format
    with open('CV_MR_constitutive.txt', 'w') as f:
        f.write(str(CV_MR_constitutive))

#!/usr/bin/env python3
"""
Purpose: Stochastic simulations constitutive model
"""
import argparse
from multiprocessing import Pool, cpu_count
import numpy as np
from scipy.stats import iqr


def gillespie_const(tmax, dt, phiR, phiM, n, R0, M0):
    """Monte-Carlo simulation of the constitutive model"""
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
        tau = np.log(1.0 / np.random.rand()) / p_t
        t = t + tau

        # Selecting reaction
        r = np.random.rand()
        cumulative = 0.0
        if r < (cumulative := cumulative + p_plus / p_t):
            p += 1
        elif r < (cumulative := cumulative + p_minus / p_t):
            p -= 1
        elif r < (cumulative := cumulative + R_plus / p_t):
            R += 1
        elif r < (cumulative := cumulative + R_minus / p_t):
            R -= 1
        elif r < (cumulative := cumulative + M_plus / p_t):
            M += 1
        else:
            M -= 1

    # If t goes beyond multiple increments of t_save, record all intermediate steps
    while t >= t_save and t_save <= tmax:
        output.append([t_save, R, M, p])
        t_save += dt

    return np.array(output)


# --------------------------------------------------


def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Stochastic simulations constitutive model",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p",
        "--parameters",
        help="csv file simulation parameters",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/constitutive_simulation_params.csv",
    )

    parser.add_argument(
        "-s",
        "--seed",
        help="set random seed for reproducibility",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/saved_seed_constitutive.txt",
    )

    parser.add_argument(
        "-w",
        "--workers",
        help="number of workers for parallel computation",
        metavar="integer",
        type=int,
        default=-1,
    )

    parser.add_argument(
        "-cv",
        "--cv_mr",
        help="output txt file with M/R CV value",
        metavar="FILE",
        type=argparse.FileType("w"),
        default="./simulations_data/CV_MR_constitutive.txt",
    )

    parser.add_argument(
        "-mr",
        "--mr_vec",
        help="output csv file with M/R simulation output",
        metavar="FILE",
        type=argparse.FileType("w"),
        default="./simulations_data/MR_constitutive_vec.csv",
    )

    return parser.parse_args()


# -------------------------------------------------------------------------------
def run_sim(_):
    """Paralelizing Monte-Carlo simulations"""
    return gillespie_const(tmax, dt, phiLambdaR, phiLambdaM, n, R0, M0)


# --------------------------------------------------------------------------------

tmax = 12
R0 = 0
M0 = 0
ts_e = 10
N = 100  # number of points recorded in the simulation
RunsDist = 2000
dt = tmax / N

# --------------------------------------------------------------------------------

if __name__ == "__main__":

    args = get_args()

    constitutive_simulation_params = args.parameters
    saved_seed_constitutive = args.seed
    workers = args.workers
    CV_MR = args.cv_mr
    MR_vec = args.mr_vec

    params = np.loadtxt(constitutive_simulation_params, delimiter=",")
    phiLambdaR, phiLambdaM, n = params

    loaded_seed = int(saved_seed_constitutive.read())
    # Set the seed
    np.random.seed(loaded_seed)

    max_workers = cpu_count()
    if not 1 <= workers <= max_workers:
        workers = max_workers

    with Pool(processes=workers) as pool:
        Simulation_constitutive = pool.map(run_sim, range(RunsDist))

    MR_vec_e = []
    # analyze simulations
    for i in range(RunsDist):
        output = Simulation_constitutive[i]
        ts = output[:, 0]
        R_constitutive = output[:, 1]
        M_constitutive = output[:, 2]

        # Find first index where ts > ts_e
        indices = np.where(ts > ts_e)[0]
        index_e = indices[0]

        if R_constitutive[index_e] > 0 and M_constitutive[index_e] > 0:
            MR_vec_e.append(M_constitutive[index_e] / R_constitutive[index_e])

    MR_vec_e = np.array(MR_vec_e)
    MR_Median = np.median(MR_vec_e)
    M_R_iqr = iqr(MR_vec_e)
    CV_MR_constitutive = M_R_iqr / MR_Median

    # Save calculated M/R value
    CV_MR.write(str(CV_MR_constitutive))
    # Save M/R values from each simulation in CSV format
    np.savetxt(MR_vec, MR_vec_e, delimiter=",")

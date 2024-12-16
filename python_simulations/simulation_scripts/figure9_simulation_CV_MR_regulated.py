#!/usr/bin/env python3
"""
Purpose: Stochastic simulations regulated model
"""
import argparse
from multiprocessing import Pool, cpu_count
import numpy as np
from scipy.stats import iqr

n = 29.915
alpha = 16.9
pC = 25
omega = 130
Kd1 = 1600
gamma = 5.1


def occupancyCR(C):
    """Calculates occupancy of P.CR promoter"""
    x = alpha / 4 * (np.sqrt(1 + 8 * C / Kd1) - 1)
    return x**2 / (1 + (1 + 1 / pC) * x**2 + omega / pC * x**4)


def occupancyM(C):
    """Calculated occupancy of P.M promoter"""
    x = alpha / 4 * (np.sqrt(1 + 8 * C / Kd1) - 1)
    return 1 / (x**2 * gamma + 1)


def gillespie_regulated(tmax, dt, params, n, R0, M0):
    """Monte-Carlo simulation of the regulated model"""
    R = R0
    p = int(round(n))
    M = M0
    phiL_lambda, phiM_lambda, beta, phiMet_lambda = params

    t = 0.0
    t_save = dt
    output = []
    output.append([t, R, M, p])

    while t < tmax:
        # Propensities
        p_plus = p
        p_minus = p
        R_plus = p * beta * phiL_lambda + p * beta * phiM_lambda * occupancyCR(R / beta)
        R_minus = R
        M_plus = p * phiMet_lambda * occupancyM(R / beta)
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
        description="Stochastic simulations regulated model",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p",
        "--parameters",
        help="csv with simulation parameter values",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/parameters_stochastic.csv",
    )

    parser.add_argument(
        "-s",
        "--seed",
        help="set random seed for reproducibility",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/saved_seed_regulated.txt",
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
        default="./simulations_data/CV_MR_regulated.txt",
    )

    parser.add_argument(
        "-mr",
        "--mr_mat",
        help="output csv file with M/R, M and R values as columns",
        metavar="FILE",
        type=argparse.FileType("w"),
        default="./simulations_data/MR_regulated_mat.csv",
    )

    return parser.parse_args()


# -------------------------------------------------------------------------------
def run_sim(_):
    """Paralelizing Monte-Carlo simulations"""
    return gillespie_regulated(tmax, dt, params, n, R0, M0)


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

    regulated_simulation_params = args.parameters
    saved_seed_regulated = args.seed
    workers = args.workers
    CV_MR = args.cv_mr
    MR_mat = args.mr_mat

    params = np.loadtxt(regulated_simulation_params, delimiter=",")
    # phiLambdaR, phiLambdaM, n = params

    loaded_seed = int(saved_seed_regulated.read())
    # Set the seed for reproducibility
    np.random.seed(loaded_seed)

    max_workers = cpu_count()
    if not 1 <= workers <= max_workers:
        workers = max_workers

    with Pool(processes=workers) as pool:
        Simulation_regulated = pool.map(run_sim, range(RunsDist))

    MR_vec_e = []
    R0_start_vec = []
    M0_start_vec = []
    # analyze simulations
    for i in range(RunsDist):
        output = Simulation_regulated[i]
        ts = output[:, 0]
        R_regulated = output[:, 1]
        M_regulated = output[:, 2]

        # Find first index where ts > ts_e
        indices = np.where(ts > ts_e)[0]
        index_e = indices[0]

        if R_regulated[index_e] > 0 and M_regulated[index_e] > 0:
            MR_vec_e.append(M_regulated[index_e] / R_regulated[index_e])
            R0_start_vec.append(R_regulated[index_e])
            M0_start_vec.append(R_regulated[index_e])

    MR_vec_e = np.array(MR_vec_e)
    R0_start_vec = np.array(R0_start_vec)
    M0_start_vec = np.array(M0_start_vec)

    MR_Median = np.median(MR_vec_e)
    M_R_iqr = iqr(MR_vec_e)
    CV_MR_regulated = M_R_iqr / MR_Median

    # Save calculated M/R value
    CV_MR.write(str(CV_MR_regulated))

    assert len(MR_vec_e) == len(R0_start_vec) == len(M0_start_vec)
    # Save M/R value as the first column, R as the second column, and M as the third column
    MR_mat_0_e = np.column_stack((MR_vec_e, R0_start_vec, M0_start_vec))
    np.savetxt(MR_mat, MR_mat_0_e, delimiter=",")

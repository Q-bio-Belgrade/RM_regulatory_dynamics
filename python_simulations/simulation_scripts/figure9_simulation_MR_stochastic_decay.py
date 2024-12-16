#!/usr/bin/env python3
"""
Purpose: Post segregation stochastic simulation
Motivation: No," Ramanujan replied, "it is a very interesting number;
it is the smallest number that can be expressed as the sum of two cubes in two distinct ways.
The Hardy-Ramanujan number: 1729 == 1**3 + 12**3 == 9**3 + 10**3
HAPPY MODELING!
"""
import argparse
from multiprocessing import Pool, cpu_count
import numpy as np
from scipy.stats import iqr


def Gillespie_decay(args):
    """Monte Carlo simulation of RM post-segregation dynamics"""
    tmax, dt, R0, M0 = args
    R = R0
    M = M0
    t = 0.0
    t_save = dt

    output = []
    output.append([t, R, M])

    while t < tmax:
        R_minus = R
        M_minus = M

        p_t = R_minus + M_minus
        if p_t == 0:
            # System static, fill up remaining times
            while t_save <= tmax:
                output.append([t_save, R, M])
                t_save += dt
            break

        # Time to next reaction
        tau = np.log(1 / np.random.rand()) / p_t

        # Select reaction
        r = np.random.rand()
        if r < (M_minus / p_t):
            # Decrease M by 1
            M -= 1
        else:
            # Decrease R by 1
            R -= 1

        t += tau

        # Record states at each dt
        while t >= t_save and t_save <= tmax:
            output.append([t_save, R, M])
            t_save += dt

    return np.array(output)

    # --------------------------------------------------


def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Post segregation stochastic simulation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--initial_state",
        help="csv with M and R initial state as columns",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/MR0_start.csv",
    )

    parser.add_argument(
        "-s",
        "--seed",
        help="set random seed for reproducibility",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/saved_seed_decay.txt",
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
        "-mr",
        "--mr_mat",
        help="output csv file with C.V M/R and median R values as columns",
        metavar="FILE",
        type=argparse.FileType("w"),
        default="./simulations_data/CVmr_R_matrix.csv",
    )

    return parser.parse_args()


if __name__ == "__main__":

    args = get_args()

    initial_state_simulation = args.initial_state
    saved_seed_decay = args.seed
    workers = args.workers
    CVmr_R_matrix = args.mr_mat

    # Read initial state from CSV
    M0_R0_start_vec = np.loadtxt(initial_state_simulation, delimiter=",")
    R0_start_vec = M0_R0_start_vec[:, 1]
    M0_start_vec = M0_R0_start_vec[:, 0]

    # *************************************************************************
    MR0_start_vec = M0_start_vec / R0_start_vec
    MR0_Median = np.median(MR0_start_vec)
    MR0_iqr = iqr(MR0_start_vec)
    CV_MR0 = MR0_iqr / MR0_Median
    CV_MR_divisions_vec = [CV_MR0]

    # ***************************************************************************
    R0_Median = np.median(R0_start_vec)
    R_divisions_vec = [R0_Median]

    # ******************************************************************************
    N = 200
    RunsDist = len(R0_start_vec)
    tmax = 8
    dt = tmax / N

    # Prepare arguments for parallel runs
    args_list = [(tmax, dt, R0_start_vec[i], M0_start_vec[i]) for i in range(RunsDist)]

    #####################################################################

    loaded_seed = int(saved_seed_decay.read())
    # Set the seed for reproducibility
    np.random.seed(loaded_seed)

    #####################################################################
    max_workers = cpu_count()
    if not 1 <= workers <= max_workers:
        workers = max_workers

    # Run simulations in parallel
    with Pool(processes=workers) as pool:
        Simulation_decay = pool.map(Gillespie_decay, args_list)

    ts_e_vec = [1, 2, 3, 4, 5, 6, 7]

    for ts_e in ts_e_vec:
        R_vec_e = []
        MR_vec_e = []
        for i in range(RunsDist):
            output = Simulation_decay[i]
            ts = output[:, 0]
            R_decay = output[:, 1]
            M_decay = output[:, 2]

            # Find first index where ts > ts_e
            indices = np.where(ts > ts_e)[0]
            index_e = indices[0]
            if R_decay[index_e] > 0 and M_decay[index_e] > 0:
                R_vec_e.append(R_decay[index_e])
                MR_vec_e.append(M_decay[index_e] / R_decay[index_e])

        R_Median = np.median(R_vec_e)
        MR_Median = np.median(MR_vec_e)
        MR_iqr = iqr(MR_vec_e)
        CV_MR = MR_iqr / MR_Median
        CV_MR_divisions_vec.append(CV_MR)
        R_divisions_vec.append(R_Median)

    ts_e_vec = [0] + ts_e_vec

    # Saving data
    CVmr_R_mat = np.column_stack([ts_e_vec, CV_MR_divisions_vec, R_divisions_vec])
    np.savetxt(CVmr_R_matrix, CVmr_R_mat, delimiter=",")

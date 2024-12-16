#!/usr/bin/env python3
"""
Purpose: R dynamics parameter inference
"""
import argparse
from multiprocessing import cpu_count
import numpy as np
from scipy.optimize import differential_evolution
from scipy.integrate import odeint

# --------------------------------------------------


# ODE function equivalent to R_loop_late
def R_loop(R, t, s, N, lambda_R, beta):
    """Calculate derivative for diff. equation solving"""
    alpha = 16.9
    p = 25
    omega = 130

    x = alpha / 4 * (np.sqrt(1 + 8 * R / beta) - 1)
    dRdt = (
        N * (s + x**2 / (1 + (1 + 1 / p) * x**2 + (omega / p) * x**4))
        - (1 / 320 + lambda_R) * R
    )
    return dRdt


# --------------------------------------------------


# integrates the ODE and then interpolates
def predicted_R(time_exp_zero, R0, params):
    """Calculate R at experimental time-points"""
    s, N, lambda_R, beta = params
    tspan = np.arange(0, max(time_exp_zero) + 2)

    R_sol = odeint(R_loop, R0, tspan, args=(N, s, lambda_R, beta))
    R_sol = R_sol.flatten()

    # Interpolate at the experimental times
    R_pred = np.interp(time_exp_zero, tspan, R_sol)
    return R_pred


# Objective function: sum of squared differences
def objective(params):
    """Sum of squared differences between experiment and predictions"""
    R_mod = predicted_R(time_exp_zero, R0, params)
    return np.sum((R_mod - R_exp) ** 2)


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="R dynamics parameter inference for Supplement Figure S1",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-e",
        "--experimental",
        help="csv file experimental data: first column time; second column R amount",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/experimental_data_R.csv",
    )

    parser.add_argument(
        "-s",
        "--seed",
        help="set random seed for reproducibility",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/saved_seed.txt",
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
        "-o",
        "--output",
        help="csv file with inferred parameters (s first parameter)",
        metavar="FILE",
        type=argparse.FileType("w"),
        default="./simulations_data/optimal_params.csv",
    )

    return parser.parse_args()


# --------------------------------------------------

args = get_args()

experimental_data = args.experimental
saved_seed = args.seed
num_workers = args.workers
output_params = args.output

# load experimental data
exp_data = np.loadtxt(experimental_data, delimiter=",", skiprows=1)
time_exp_zero = exp_data[:, 0]
R_exp = exp_data[:, 1]
R0 = R_exp[0]

# set random seed for reproducibility
seed_value = int(saved_seed.read())
np.random.seed(seed_value)

max_workers = cpu_count()
if not 1 <= num_workers <= max_workers:
    num_workers = max_workers

lb = [0.01, 200, 1 / 5000, 1 / 50000]
ub = [0.3, 400, 1, 1]
bounds = list(zip(lb, ub))

# Run the GA-like optimizer
result = differential_evolution(
    objective,
    bounds,
    maxiter=4000,
    tol=1e-10,
    updating="deferred",
    workers=max_workers,  # parallelization
    seed=seed_value,
)

optimalParams = result.x
fval = result.fun

np.savetxt(output_params, optimalParams, delimiter=",")

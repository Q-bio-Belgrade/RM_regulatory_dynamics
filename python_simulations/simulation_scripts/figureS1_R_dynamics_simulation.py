#!/usr/bin/env python3
"""
Purpose: R dynamics simulation
"""

import argparse
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="R dynamics simulation, Supplement Figure S1",
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
        "-t",
        "--time_first",
        help="initial time for plotting R data",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/initial_time.txt",
    )

    parser.add_argument(
        "-p",
        "--parameters",
        help="csv file with inferred simulation parameters",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/optimal_params.csv",
    )

    return parser.parse_args()


# --------------------------------------------------
def main():

    args = get_args()

    experimental_data = args.experimental
    t_first = args.time_first
    loaded_params = args.parameters

    # load experimental data
    exp_data = np.loadtxt(experimental_data, delimiter=",", skiprows=1)
    time_exp_zero = exp_data[:, 0]
    R_exp = exp_data[:, 1]
    R0 = R_exp[0]

    params = np.loadtxt(loaded_params, delimiter=",")

    # The function calculating dRdt
    def R_loop_late(t, R, s, N, lambda_R, beta):
        """Calculate derivative for diff. equation solving"""
        alpha = 16.9
        p = 25
        omega = 130

        x = alpha / 4 * (np.sqrt(1 + 8 * R / beta) - 1)

        numerator = x**2
        denominator = 1 + (1 + 1 / p) * x**2 + (omega / p) * x**4
        dRdt = N * (s + numerator / denominator) - (1 / 320 + lambda_R) * R
        return dRdt

    # The function calculating simulated time and R
    def R_simulation(time_exp_zero, R0, params):
        """Simulate R dynamics"""
        s, N, lambda_R, beta = params

        t_max = np.max(time_exp_zero) + 1
        t_eval = np.arange(0, t_max + 1, 1)

        # Use solve_ivp to solve the ODE
        sol = solve_ivp(
            lambda t, R: R_loop_late(t, R, s, N, lambda_R, beta),
            [0, t_max],
            [R0],
            t_eval=t_eval,
        )

        T_zero = sol.t
        R = sol.y[0]
        return T_zero, R

    # Perform the simulation
    T_zero, R = R_simulation(time_exp_zero, R0, params)

    time_first = np.loadtxt(t_first)

    time_exp = time_exp_zero + time_first
    T = T_zero + time_first

    # Plot the results
    plt.plot(T, R, "-", label="Simulation")
    plt.plot(time_exp, R_exp, "*", label="Experimental")

    plt.xlabel("time (min)", fontsize=15)
    plt.ylabel("R amount (relative units)", fontsize=15)

    ax = plt.gca()
    ax.tick_params(direction="out", length=5, width=1, labelsize=14)

    plt.legend()
    plt.grid(True)
    plt.savefig("./simulation_results/FigureS1.png", dpi=300)


# --------------------------------------------------
if __name__ == "__main__":
    main()

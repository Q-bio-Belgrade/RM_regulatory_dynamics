#! /usr/bin/env python3
"""
Purpose: Panel plot of stochastic simulation results
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Panel plot of stochastic simulation results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-r",
        "--regulated_output",
        help="csv file with regulated model stochastic simulation output",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/MR_regulated_mat.csv",
    )

    parser.add_argument(
        "-c",
        "--constitutive_output",
        help="csv file with constitutive model stochastic simulation output",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/MR_constitutive_vec.csv",
    )

    parser.add_argument(
        "-f",
        "--cv_regulated",
        help="calculated M/R C.V. for regulated model",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/CV_MR_regulated.txt",
    )

    parser.add_argument(
        "-b",
        "--cv_constitutive",
        help="calculated M/R C.V. for constitutive model",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/CV_MR_constitutive.txt",
    )

    parser.add_argument(
        "-d",
        "--decay_output",
        help="csv file with post-segregation killing simulation output",
        metavar="FILE",
        type=argparse.FileType("rt"),
        default="./simulations_data/CVmr_R_matrix.csv",
    )

    parser.add_argument(
        "-o",
        "--figure_panel",
        help="Saved panel for stochastic simulation results(Fig. 9)",
        metavar="FILE",
        type=str,
        default="./simulation_results/Figure9.png",
    )

    return parser.parse_args()


# --------------------------------------------------
def main():

    args = get_args()

    regulated_output = args.regulated_output
    constitutive_output = args.constitutive_output
    CV_regulated = args.cv_regulated
    CV_constitutive = args.cv_constitutive
    decay_output = args.decay_output
    figure = args.figure_panel

    # Read MR_regulated_mat and compute CV
    MR_regulated_mat = pd.read_csv(regulated_output, header=None).values
    MR_vec_e = MR_regulated_mat[:, 0]

    # Initialize figure and a grid layout
    fig = plt.figure(figsize=(10, 6))
    gs = fig.add_gridspec(2, 5)

    # First subplot
    # Spanning 2 rows and 2 columns: columns 0-1
    ax1 = fig.add_subplot(gs[:, 0:2])

    # Plot histogram for regulated (blue)
    ax1.hist(MR_vec_e, facecolor="b", alpha=0.7, edgecolor="black", bins=40)
    ax1.set_xlim([0, 0.9])
    ax1.set_xlabel("M/R")
    ax1.set_ylabel("simulated distributions")
    ax1.hold = True

    # Read MR_constitutive_vec and compute CV
    MR_constitutive_vec = pd.read_csv(constitutive_output, header=None).values.flatten()

    # Overplot histogram for constitutive (red)
    ax1.hist(MR_constitutive_vec, facecolor="r", alpha=0.7, edgecolor="black", bins=40)

    # Add legend
    labels = ["regulated", "constitutive"]
    ax1.legend(labels, loc="upper left")

    # Second subplot
    # Spanning 2 rows and 1 column: column 2
    ax2 = fig.add_subplot(gs[:, 2])

    x = ["constitutive", "regulated", "experimental"]
    CV_MR_experimental = 0.3918

    # Read CV_MR values from files
    CV_MR_constitutive = np.loadtxt(CV_constitutive)
    CV_MR_regulated = np.loadtxt(CV_regulated)

    bar_values = [CV_MR_constitutive, CV_MR_regulated, CV_MR_experimental]

    b = ax2.bar(x, bar_values, color="gray")

    ax2.set_xticklabels(x, rotation=45, ha="right")

    # Set individual colors for each bar
    b[0].set_color([1, 0, 0])  # red
    b[1].set_color([0, 0, 1])  # blue
    b[2].set_color([0, 1, 0])  # green

    ax2.set_ylabel("coefficient of variation")

    # Third subplot
    # Spanning 2 rows and 2 columns: columns 3-5
    ax3 = fig.add_subplot(gs[:, 3:5])

    # Read the CSV file into a DataFrame
    CVmr_R_mat = pd.read_csv(decay_output, header=None)

    # Extract columns as numpy arrays
    ts_e_vec = CVmr_R_mat.iloc[:, 0].to_numpy()
    CV_MR_vec = CVmr_R_mat.iloc[:, 1].to_numpy()
    R_Median_vec = CVmr_R_mat.iloc[:, 2].to_numpy()

    # Create dual y-axis
    ax3_right = ax3.twinx()

    # Right y-axis: CV M/R (red)
    ax3_right.plot(ts_e_vec, CV_MR_vec, "-r", linewidth=2, label="CV M/R")
    ax3_right.set_ylabel("CV M/R")
    ax3_right.set_ylim([0, 1.2])

    # Left y-axis: M/R median (blue dashed)
    ax3.plot(ts_e_vec, R_Median_vec, "b--", linewidth=2, label="M/R median")
    ax3.set_ylabel("R molecules")
    ax3.set_xlabel("cell divisions")

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(figure, dpi=300)
    plt.show()


# --------------------------------------------------
if __name__ == "__main__":
    main()

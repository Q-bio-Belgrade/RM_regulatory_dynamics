import time
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import iqr

def Gillespie_decay(args):
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


if __name__ == '__main__':
    # Read initial state from CSV
    M0_R0_start_vec = np.loadtxt("MR0_start.csv", delimiter=',')
    R0_start_vec = M0_R0_start_vec[:, 1]
    M0_start_vec = M0_R0_start_vec[:, 0]

    #*************************************************************************
    MR0_start_vec = M0_start_vec / R0_start_vec
    MR0_Median = np.median(MR0_start_vec)
    MR0_iqr = iqr(MR0_start_vec)
    CV_MR0 = MR0_iqr / MR0_Median
    CV_MR_divisions_vec = [CV_MR0]

    #***************************************************************************
    R0_Median = np.median(R0_start_vec)
    R_divisions_vec = [R0_Median]

    #******************************************************************************
    N = 200
    RunsDist = len(R0_start_vec)
    tmax = 8
    dt = tmax / N

    # Prepare arguments for parallel runs
    args_list = [(tmax, dt, R0_start_vec[i], M0_start_vec[i]) for i in range(RunsDist)]

    #####################################################################

    #Generate a seed based on the current time
    seed_value = int(time.time())
    np.random.seed(seed_value)

    # Save this seed to a file
    with open("saved_seed_decay.txt", "w") as f:
        f.write(str(seed_value))

    #with open("saved_seed_constitutive.txt", "r") as f:
        #loaded_seed = int(f.read())

    # Set the seed
    #np.random.seed(loaded_seed)

    #####################################################################


    # Run simulations in parallel
    with Pool(processes=8) as pool:
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

    # Plotting
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax2.plot(ts_e_vec, CV_MR_divisions_vec, '-r', linewidth=2, label="M/R C.V.")
    ax2.set_ylabel('C.V. M/R')
    ax2.set_ylim([0, 1.2])

    ax1.plot(ts_e_vec, R_divisions_vec, 'b--', linewidth=2, label="R median")
    #ax1.set_ylim([0, 1.2e4])
    ax1.set_ylabel('R (molecules)')
    ax1.set_xlabel('Time (cell divisions)')

    # If you want to display legend
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc='upper right')

    plt.tight_layout()
    plt.savefig("plot_output.png", dpi=300)
    plt.show()

    # Saving data
    CVmr_R_matrix = np.column_stack([ts_e_vec, CV_MR_divisions_vec, R_divisions_vec])
    np.savetxt("CVmr_R_matrix.csv", CVmr_R_matrix, delimiter=',')

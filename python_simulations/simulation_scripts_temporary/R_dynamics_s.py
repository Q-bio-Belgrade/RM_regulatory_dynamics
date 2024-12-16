#! /usr/bin/env python3

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Simulation parameters in the format [ N,s,lambda_R,beta ]
optimalParams = [274.7185, 0.2099, 0.0002, 314.6847]
N     = optimalParams[0]
s     = optimalParams[1]
lambda_R = optimalParams[2]
beta  = optimalParams[3]

# load experimental data
exp_data = np.loadtxt("./simulations_data/experimental_data_R.csv",delimiter=",",skiprows=1)
time_exp_zero = exp_data[:,0]
R_exp = exp_data[:,1]
R0 = R_exp[0]

# The function calculating dRdt
def R_loop_late(t, R, N, s, lambda_R, beta):
    alpha = 16.9
    p = 25
    omega = 130

    x = alpha/4 * ( np.sqrt(1 + 8*R/beta) - 1 )

    # Compute dRdt
    numerator = x**2
    denominator = 1 + (1 + 1/p)*x**2 + (omega/p)*x**4
    dRdt = N*(s + numerator/denominator) - (1/320 + lambda_R)*R
    return dRdt

# The function calculating simulated time and R
def R_simulation(time_exp_zero, R0, params):
    N, s, lambda_R, beta = params

    t_max = np.max(time_exp_zero) + 1
    t_eval = np.arange(0, t_max+1, 1)  # integer steps like MATLAB

    # Use solve_ivp to solve the ODE
    sol = solve_ivp(lambda t, R: R_loop_late(t, R, N, s, lambda_R, beta),
                    [0, t_max], [R0], t_eval=t_eval)

    T_late = sol.t
    R_late = sol.y[0]
    return T_late, R_late

# Perform the simulation
T_zero, R = R_simulation(time_exp_zero, R0, optimalParams)

time_first = 408.5
time_exp = time_exp_zero + time_first
T = T_zero + time_first

# Plot the results
plt.plot(T, R, '-', label="Simulation")
plt.plot(time_exp, R_exp, '*', label="Experimental")

plt.xlabel("time (min)", fontsize=15)
plt.ylabel("R amount (relative units)", fontsize=15)

ax = plt.gca()
ax.tick_params(direction='out', length=5, width=1, labelsize=14)

plt.legend()
plt.grid(True)
plt.show()

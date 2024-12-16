import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.optimize import differential_evolution
from scipy.integrate import odeint


seed_value = int(time.time())
np.random.seed(seed_value)

#with open("saved_seed.txt", "w") as f:
#    f.write(str(seed_value))

#with open("saved_seed.txt", "r") as f:
#    seed_value = int(f.read())

np.random.seed(seed_value)

# Experimental data (from the MATLAB code)
time_exp_zero = np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330])
R_exp = np.array([3402, 4182, 5536, 6329, 6936, 8061, 9334, 9991, 10904, 13224, 11620, 13989])

time_first = 408.5
R0 = R_exp[0]

# Parameter bounds: [N, s, lambda_R, beta]
lb = [200, 0.01, 1/5000, 10000]
ub = [400, 0.3, 1, 100000]
bounds = list(zip(lb, ub))

# Define the ODE function equivalent to R_loop_late
def R_loop_late(R, t, N, s, lambda_R, beta):
    alpha = 16.9
    p = 25
    omega = 130

    # Compute x as in MATLAB code:
    x = alpha/4 * ( np.sqrt(1 + 8*R/beta) - 1 )
    dRdt = N*( s + x**2/(1+(1+1/p)*x**2 + (omega/p)*x**4 ) ) - (1/320 + lambda_R)*R
    return dRdt

# predicted_R function: integrates the ODE and then interpolates
def predicted_R(time_exp_zero, R0, params):
    N, s, lambda_R, beta = params
    tspan_late = np.arange(0, max(time_exp_zero)+2)  # +1 as in MATLAB, +1 again to ensure coverage

    R_sol = odeint(R_loop_late, R0, tspan_late, args=(N, s, lambda_R, beta))
    R_sol = R_sol.flatten()  # odeint returns Nx1 array, flatten to 1D

    # Interpolate at the experimental times
    R_pred = np.interp(time_exp_zero, tspan_late, R_sol)
    return R_pred

# R_simulation function: just returns the integrated solution over the same tspan
def R_simulation(time_exp_zero, R0, params):
    N, s, lambda_R, beta = params
    tspan_late = np.arange(0, max(time_exp_zero)+2)  # same logic as above

    R_sol = odeint(R_loop_late, R0, tspan_late, args=(N, s, lambda_R, beta))
    R_sol = R_sol.flatten()

    return tspan_late, R_sol

# Objective function: sum of squared differences
def objective(params):
    R_mod = predicted_R(time_exp_zero, R0, params)
    return np.sum((R_mod - R_exp)**2)

# Run the GA-like optimizer
result = differential_evolution(
    objective,
    bounds,
    maxiter=40000,
    tol=1e-10,
    updating='deferred',
    workers=-1,    # parallelization if available
    seed=seed_value
)

optimalParams = result.x
fval = result.fun

print(f"N: {optimalParams[0]}")
print(f"s: {optimalParams[1]}")
print(f"lambdaR: {optimalParams[2]}")
print(f"beta: {optimalParams[3]}")

T, R = R_simulation(time_exp_zero, R0, optimalParams)
time_exp = time_exp_zero + time_first
T_shift = T + time_first

plt.figure()
plt.plot(time_exp, R_exp, '*', label='Data')
plt.plot(T_shift, R, '-', label='Model')
plt.xlabel('Time')
plt.ylabel('R')
plt.legend()
plt.show()

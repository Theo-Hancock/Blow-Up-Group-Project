#%% intersection plot 1
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 3, 300)

lambdas = [1, np.e, 5]
colors = ['blue', 'red', 'green']
labels = [r'$y = \lambda \beta T, \lambda = 1$', r'$\lambda = e$', r'$\lambda = 6$']

plt.figure(figsize=(9, 6))
plt.plot(x, np.exp(x), 'k-', linewidth=2, label=r'$y = e^{\beta T}$')

for lam, color, label in zip(lambdas, colors, labels):
    plt.plot(x, lam * x, color=color, linestyle='-', linewidth=1.5, label=label)

plt.xlabel(r'$\beta T$', fontsize=12)
plt.ylabel(r'$y$', fontsize=12)
plt.legend(loc='upper left', fontsize=12)
# plt.grid(False, alpha=0.3)
plt.xlim(0, 3)
plt.ylim(0, 20)
plt.show()


#%% bifurcation diagram 1
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 5, 500)

plt.figure(figsize=(10, 5))

plt.scatter(np.e, 1, color='red', label='(e, 1)')
# 绘制不稳定区域（虚线）
x_s = x[x < 1]
plt.plot(np.exp(x_s) / x_s, x_s, 'k-', linewidth=2, label='stable')

x_un = x[x > 1]
plt.plot(np.exp(x_un) / x_un, x_un, 'k--', linewidth=2, label='unstable')

plt.xlabel(r'$\lambda$', fontsize=12)
plt.ylabel(r'$\beta T$', fontsize=12)
plt.xlim(0, 20)

plt.legend(fontsize=12)
plt.show()



#%% bifurcation diagram 2
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def eqn(T_0, T_f, beta, c, a, D_T):
    return T_0 - 2 / beta * np.log(np.cosh(np.sqrt(c * a * beta / (2 * D_T)) * np.exp(beta * T_0 / 2))) - T_f

def f(T_f, beta, c, a, D_T):
    # Solve for T_0 using fsolve
    T_0_initial_guess = T_f  # Initial guess for T_0
    T_0_solution = fsolve(eqn, T_0_initial_guess, args=(T_f, beta, c, a, D_T))

    return T_0_solution[0]


# Example usage
T_f, beta, c, a, D_T = 1.0, 0.5, 2.0, 1.0, 0.1
T_0 = f(T_f, beta, c, a, D_T)
print("Solved T_0:", T_0)

x = np.linspace(0, 300, 600)

plt.figure(figsize=(10, 5))

def f(T_f, beta, c, a, D_T):
    T_0 - 2 / beta * np.log(np.cosh(np.sqrt(c * a * beta / D_T / 2) * np.exp(beta * T_0 / 2))) = T_f
    return T_0

# plt.scatter(np.e, 1, color='red', label='(e, 1)')
# 绘制不稳定区域（虚线）
x_s = x[x < 1]
plt.plot(np.exp(x_s) / x_s, x_s, 'k-', linewidth=2, label='stable')

x_un = x[x > 1]
plt.plot(np.exp(x_un) / x_un, x_un, 'k--', linewidth=2, label='unstable')

plt.xlabel(r'$\lambda$', fontsize=12)
plt.ylabel(r'$\beta T$', fontsize=12)
plt.xlim(0, 20)

plt.legend(fontsize=12)
plt.show()

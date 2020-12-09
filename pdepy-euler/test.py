import numpy as np
import matplotlib.pyplot as plt

'''
    EXAMPLE
Centered Euler method for u' = ru; u(0) = u_0 on interval [0, 50]
'''


u_0 = 10
r = 0.2
dx = 1
N = 50

x = np.linspace(0, (N + 1) * dx, N + 2)
u = np.zeros(N + 2)

# initial condition
u[0] = u_0
# forward euler approximation for u(0 + dx)
u[1] = (r * u[0] * dx) + u[0]

# centered euler approximation for u(x)
for n in range(2, N + 2):
    u[n] = (r * u[n - 1] * 2 * dx) + u[n - 2]

plt.plot(x, u, 'b-')
plt.plot(x, u_0 * np.exp(r * x), 'r-')
plt.legend(['numerical', 'exact'], loc = 'upper left')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.savefig('test centered euler 1')
plt.show()

plt.close()


'''
IVP N''(t) = r * N(t) ; N(0) = N_0 and N'(0) = N_1
Centered Euler method
'''

'''
N_0 = 1000
N_1 = 1
r = 3
dt = 0.001
N_t = 50000

t = np.linspace(0, (N_t + 1) * dt, N_t + 2)
N = np.zeros(N_t + 2)
N_pm = np.zeros(N_t + 2)
N[0] = N_0
N_pm[0] = N_1
# forward Euler to get N(0 + dt)
N[1] = N_pm[0] * dt + N[0]
# centered Euler approx'n for N(t)
for i in range(2, N_t + 2):
    N[i] = r * dt * dt * N[i - 1] - N[i - 2] + 2 * N[i - 1]

plt.plot(t, N, 'b-')
plt.xlabel('t')
plt.ylabel('N(t)')
plt.savefig('centereuler2')
plt.show()
'''

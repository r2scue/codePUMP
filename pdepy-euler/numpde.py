import numpy as np
from scipy import integrate
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


### EX: general 1-dim heat eqn
#### forward-time centered-space scheme
def heateqn(k, g, I, A, B, dt, T, dx, L):
    """
    IVP for u(x, t):
        u_t = k * u_xx + g(x, t); x[i] on [0, L]
        u(x, 0) = I(x); u(0, t) = A(t); u(L, t) = B(t)

    PARAMETERS:
        k = const coefficient
        g = inhomogeneous part (fn of x, t in general)
        I = initial cond'n
        A = boundary cond'n at x = 0
        B = boundary cond'n at x = L
        dt = time step size
        T = end time
        dx = x step size
        L = upper boundary of x

    RETURNS numerical approximation for u(x[i], t[j]), U, (N_x + 1) by (N_t + 1)
    array
    """
    # initialize data
    N_x = int(round(float(L / dx))) # number of spatial steps (x)
    x = np.linspace(0, L, N_x + 1)
    N_t = int(round(float(T / dt))) # number of time steps (t)
    t = np.linspace(0, T, N_t + 1)

    U = np.zeros((N_x + 1, N_t + 1))
    U[:, 0] = I(x)
    U[0, :] = A(t)
    U[N_x, :] = B(t)

    s = k * dt / dx**2
    for j in range(0, N_t):
        for i in range(1, N_x - 1):
            U[i, j + 1] = s * (U[i + 1, j] - 2 * U[i, j] + U[i - 1, j]) + \
                            g(x[i], t[j]) * dt + U[i, j]

    return U, t


### EX: general 1-dim wave eqn for vibrating string
#### centered-time centered-space scheme
def vibstr(c, h, I, G, A, B, dt, T, dx, L):
    """
    IVP for u(x, t):
        u_tt = (c^2)u_xx + h(x, t); x[i] on [0, L]
        u(x, 0) = I(x); u_t(x, 0) = G(x); u(0, t) = A(t); u(L, t) = B(t)

    PARAMETERS:
        c = wave speed: const, > 0
        h = inhomogeneous part (source/sink): fn of x, t
        I = initial condition on u: fn of x
        G = initial condition on u_t: fn of x
        A = boundary condition at x = 0: fn of t
        B = boundary condition at x = L: fn of t
        dt = time step size
        T = total time
        dx = spatial step size
        L = upper boundary on x

    RETURNS:
        U: numerical sol'n for vertical displacement w.r.t. x, t
        X: 2d expansion of 1d array x
        T: 2d expansion of 1d array t
    """
    N_t = int(round(float(T / dt)))
    t = np.linspace(0, T, N_t + 1)
    N_x = int(round(float(L / dx)))
    x = np.linspace(0, L, N_x + 1)

    U = np.zeros((N_x + 1, N_t + 1))
    U[:, 0] = I(x)
    U[0, :] = A(t)
    U[N_x, :] = B(t)

    ## calculated values
    beta = (c**2) * (dt**2) /(dx**2)

    u_jeqminus1 = np.zeros(N_x + 1)
    u_jeqminus1[0] = A(-1 * dt)
    u_jeqminus1[N_x] = B(-1 * dt)

    ## populate U
    ### col j = -1, 1:
    for i in range(1, N_x):
        u_jeqminus1[i] = (beta * (I(x[i + 1]) - 2 * I(x[i]) + I(x[i - 1])) +
                            (dt**2) * h(x[i], 0) + 2 * I(x[i]) - 2 * dt * G(x[i])) / 2

        U[i, 1] = beta * (U[i + 1, 0] - 2 * U[i, 0] + U[i - 1, 0]) + \
                    (dt**2) * h(x[i], 0) + 2 * U[i, 0] - u_jeqminus1[i]

    ### cols j = 2 to N_t - 1, rows i = 1 to N_x:
    for jtr in range(1, N_t):
        for itr in range(1, N_x):
            U[itr, jtr + 1] = beta * (U[itr + 1, jtr] - 2 * U[itr, jtr] +
                                U[itr - 1, jtr]) + \
                                (dt**2) * h(x[itr], t[jtr]) + \
                                2 * U[itr, jtr] - U[itr, jtr - 1]
            print(itr, jtr)

    ## create mesh grid x, t
    X = np.transpose( np.tile(x, (N_t + 1, 1)) ) 
    T = np.tile(t, (N_x + 1, 1))

    return U, X, T


### EX: vibstr demo, vibrating string IVP on large interval
def ex1(c, dt, T, dx, L):
    """
    Numerical solution for IVP: homogeneous vibrating string, initially single 
    sin pulse at rest, Dirichlet BCs

    PARAMETERS:
        c = wave speed
        dt = time step size
        T = total time
        dx = x step size
        L = upper boundary on x (>= 16 suggested)
    """

    def I(x):
        return (np.heaviside(x-8*np.pi, 1) - np.heaviside(x-9*np.pi, 1)) * np.sin(x)

    U, x, t = vibstr(
        c = c,
        h = lambda x, t: 0,
        I = I,
        G = lambda x: 0,
        A = lambda t: 0,
        B = lambda t: 0,
        dt = dt,
        T = T,
        dx = dx,
        L = L
    )

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    #ax.contour3D(x, t, U, 50, cmap='binary')
    ax.plot_wireframe(x, t, U)
    #ax.plot_surface(x, t, U, rstride=1, cstride=1,
    #                cmap='viridis', edgecolor='none')

    ax.set_zlim3d(bottom=-1.5, top=1.5)

    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u')
    ax.set_title('ex1 CTCS: u(x, t)')

    ax.view_init(30, -120)

    plt.show()


def ex0_fs():
    c = 1
    dt = 0.01
    dx = 0.1
    T = 5
    L = 10

    N_t = int(round(float(T / dt))) # index of t = T
    t = np.linspace(0, T, N_t + 1)
    N_x = int(round(float(L / dx))) # index of x = L
    x = np.linspace(0, L, N_x + 1)

    U_fs = np.zeros((N_x + 1, N_t + 1))

    for i in range(N_x + 1):
        for j in range(N_t + 1):
            n = 1
            while n <= 50:
                U_fs[i, j] = U_fs[i, j] + \
                                (2 * L * np.sin(n * (np.pi**2) / L) * 
                                np.cos(n * np.pi * c * t[j] / L) * 
                                np.sin(n * np.pi * x[i] / L)) / \
                                (L**2 - (n * np.pi)**2)
                n += 1

    X = np.transpose( np.tile(x, (N_t + 1, 1)) ) 
    T = np.tile(t, (N_x + 1, 1))

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    #ax.contour3D(X, T, U_e, 50, cmap='binary')
    ax.plot_wireframe(X, T, U_fs, color='red')
    #ax.plot_surface(X, T, U_e, rstride=1, cstride=1,
    #                cmap='viridis', edgecolor='none')

    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u')
    ax.set_title('ex0 FS: u(x, t)')

    ax.view_init(10, -150)

    plt.show()


def ex0_c():
    c = 1
    dt = 0.1
    dx = 0.01
    T = 5
    L = 10

    U, x, t = vibstr(
        c = c,
        h = lambda x, t: 0,
        I = lambda x: (np.heaviside(x, 1) - np.heaviside(x - np.pi, 1)) * np.sin(x),
        G = lambda x: 0,
        A = lambda t: 0,
        B = lambda t: 0,
        dt = dt,
        T = T,
        dx = dx,
        L = L
    )

    N_t = int(round(float(T / dt))) # index of t = T
    t1 = np.linspace(0, T, N_t + 1)
    N_x = int(round(float(L / dx))) # index of x = L
    x1 = np.linspace(0, L, N_x + 1)

    U_fs = np.zeros((N_x + 1, N_t + 1))

    for i in range(N_x + 1):
        for j in range(N_t + 1):
            n = 1
            while n <= 50:
                U_fs[i, j] = U_fs[i, j] + \
                                (2 * L * np.sin(n * (np.pi**2) / L) * 
                                np.cos(n * np.pi * c * t1[j] / L) * 
                                np.sin(n * np.pi * x1[i] / L)) / \
                                (L**2 - (n * np.pi)**2)
                n += 1

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_wireframe(x + 0.01, t, U + 0.01)
    ax.plot_wireframe(x, t, U_fs, color='red')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u')
    ax.set_title('ex0 FS (red) and CTCS (blue): u(x, t)')
    ax.view_init(10, -150)
    plt.show()


if __name__ == '__main__':
    #ex0_c() # vibstr test against FS approximation
    ex1(1, 0.01, 10, 0.01, 50) # infinite vibstr sim


























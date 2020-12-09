import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpde import vibstr


def ex1_mv():
    c = 1
    dt = 0.01
    dx = 0.01
    T = 10
    L = 30

    def I(x):
        return (np.heaviside(x - 4 * np.pi, 1) - np.heaviside(x - 5 * np.pi, 1)) * np.sin(x)

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

    x_ = x[:, 0]
    t_ = t[0, :]

    fig = plt.figure()
    ax = plt.axes(xlim = (0, L), ylim = (-1.5, 1.5))
    line, = ax.plot([], [], lw = 2)

    def init():
        line.set_data(x[:, 0], U[:, 0])
        return line,


    def animate(i):
        i_ = int(i)
        line.set_data(x[:, i_], U[:, i_])
        line.set_color('red')
        return line, 

    anim = animation.FuncAnimation(fig, animate, init_func = init, frames = len(t_), 
            interval = 1)
    plt.show()


if __name__ == '__main__':
    ex1_mv()




















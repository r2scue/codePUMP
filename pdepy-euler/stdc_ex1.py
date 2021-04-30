import numpy as np 
#import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from mpl_toolkits import mplot3d


# const
pi = np.pi


# IC1
def f(x):
	return (np.heaviside(x - 3.45,1) - np.heaviside(x - 3.55,1)) * \
			np.sin(10*pi*(x - 0.05))


# IC2
def g(x):
	return 0


T = 3 # total time
L = 7 # total space (x direction)


# coeffs
a1 = 0.6
a2 = 0.2
b1 = 3
b2 = 1
# modulate (alpha, beta) between (a1, b1) and (a2, b2) every L/10 space step 
# and every T/10 space step

# c1=c2=c 
# time diff time mod / time diff space interface = 1/c

def alpha(x, t, m = 18, n = 210):
	xn = int(round(n * x / L))
	tm = int(round(m * t / T))

	if xn%2 == 0:
		if tm%2 == 0:
			return a1
		else:
			return a2
	else:
		if tm%2 == 0:
			return a2
		else :
			return a1


def beta(x, t, m = 18, n = 210):
	xn = int(round(n * x / L))
	tm = int(round(m * t / T))

	if xn%2 == 0:
		if tm%2 == 0:
			return b1
		else:
			return b2
	else:
		if tm%2 == 0:
			return b2
		else:
			return b1


# vars
dt_ = 0.005
dx_ = 0.005
N_t = int(round(T/dt_, 10))
N_x = int(round(L/dx_, 10))

# data
t = np.linspace(0, T, N_t+1)
x = np.linspace(0, L, N_x+1)
dt = round(t[1] - t[0], 10)
dx = round(x[1] - x[0], 10)
u = np.zeros((N_x+1, N_t+1)) # u at i=0, N_x+1 is 0
u_tm1 = np.zeros(N_x+1) # u_tm1 at i = 0, N_x+1 is 0


# IC1: compute u at j=0
u[:,0] = f(x)

# IC2: compute u_tm1 (u at j=-1) and u at j=1
for i in range(1, N_x):
	u_tm1[i] = round(
				((dt**2) / ((beta(i*dx,0.5)+beta(i*dx,-0.5))*(dx**2))) * 
				(alpha((i+0.5)*dx,0) * (f(x[i+1]) - f(x[i])) - 
				alpha((i-0.5)*dx,0) * (f(x[i]) - f(x[i-1]))) + f(i*dx) - 
				2*dt*g(i*dx)*beta(i*dx,0.5)/(beta(i*dx,0.5)+beta(i*dx,-0.5))
				, 10)

	u[i,1] = round(
			((dt**2) / (beta(i*dx,0.5)*(dx**2))) * 
			(alpha((i+0.5)*dx,0) * (f(x[i+1]) - f(x[i])) - 
			alpha((i-0.5)*dx,0) * (f(x[i]) - f(x[i-1]))) + 
			(beta(i*dx,-0.5)/beta(i*dx,0.5))*(f(x[i]) - u_tm1[i]) + 
			f(x[i])
			, 10)

# compute u at j=2 to N_t, at i=1 to N_x-1
for j in range(1, N_t):
	for i in range(1, N_x):
		u[i,j+1] = round(
					((dt**2) / (beta(i*dx,(j+0.5)*dt)*(dx**2))) * 
					(alpha((i+0.5)*dx,j*dt) * (u[i+1,j] - u[i,j]) - 
					alpha((i-0.5)*dx,j*dt) * (u[i,j] - u[i-1,j])) + 
					(beta(i*dx,(j-0.5)*dt)/beta(i*dx,(j+0.5)*dt)) * 
					(u[i,j] - u[i,j-1]) + u[i,j]
					, 10)

		print(i, j)

#		if u[i,j] != 0:
#			print(x[i], t[j], u[i,j])


fig = plt.figure()


#Writer = animation.writers['ffmpeg'] 
#writer = Writer(fps=40, metadata=dict(artist='Me'), bitrate=1800) 



ax = plt.axes(xlim = (1.25, 5.75), ylim=(-5, 5))
line, = ax.plot([], [], lw = 2)
time = ax.text(5, 0.5, s='')


def init():
	line.set_data(x, u[:, 0])
	time.set_text('initial')
	return line,


def animate(i):
	i_ = int(i)
	line.set_data(x, u[:, i_])
	time.set_text('time:'+str(t[i_]))
	line.set_color('blue')
	return line, 

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(t), 
        interval=1)

'''
ax = plt.axes(projection='3d', xlim=(0, T), ylim=(2,5), zlim=(-0.5,1))
time = ax.text(T/2, 5, 1.25, s='')
t_ = t.tolist()
x_ = x.tolist()
u_ = []
j = int(0)

def init():
	ax.scatter(t_[0], x_, u[:,0], c=u[:,0])
	time.set_text('initial')

def animate(k):
	global j
	global u_
	j = j+14

	if j >= N_t:
		j = 0
		anim.event_source.stop()
		return

	print(j)

	u_[:] = []
	for i in range(N_x+1):
		u_.append(u[i, j])

	ax.scatter(t_[j], x_, u_, c=u_, s=2, alpha=0.1)
	time.set_text('time:'+str(t[j]))


ax.set_title('STDC ex. 1 with time as space variable')
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('u')
ax.view_init(35, 180)


anim = animation.FuncAnimation(fig, animate, frames=len(t), 
		interval=1)

'''

#anim.save('st_blowup.mp4', writer=writer)

plt.show()





















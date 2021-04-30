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
	return (np.heaviside(x - 3.4,1) - np.heaviside(x - 3.5,1)) * \
			np.sin(10*pi*x)


# IC2
def g(x):
	return 0


# coeffs
def alpha(x):
	xm = np.abs(x - 3.45)

	if xm < 0.3:
		return 0.2
	elif xm >= 0.3 and xm < 0.5:
		return 0.4
	elif xm >= 0.5 and xm < 0.7:
		return 0.6
	else:
		return 1
	'''
	if xm < 0.5:
		return 0.4
	return 0.6
	'''


def beta(x):
	xm = np.abs(x - 3.45)

	if xm > 0.3 or xm < 0.5:
		return 2
	return 4
	'''
	if xm < 0.5:
		return 3
	return 6
	'''


# vars
T = 2.5
L = 7
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
				((dt**2) / (2*beta(i*dx)*(dx**2))) * 
				(alpha((i+0.5)*dx) * (f(x[i+1]) - f(x[i])) - 
				alpha((i-0.5)*dx) * (f(x[i]) - f(x[i-1]))) + 
				f(i*dx) - dt*g(i*dx)
				, 10)

	u[i,1] = round(
			((dt**2) / (beta(i*dx)*(dx**2))) * 
			(alpha((i+0.5)*dx) * (f(x[i+1]) - f(x[i])) - 
			alpha((i-0.5)*dx) * (f(x[i]) - f(x[i-1]))) + 
			2*f(x[i]) - u_tm1[i]
			, 10)

# compute u at j=2 to N_t, at i=1 to N_x-1
for j in range(1, N_t):
	for i in range(1, N_x):
		u[i,j+1] = round(
					((dt**2) / (beta(i*dx)*(dx**2))) * 
					(alpha((i+0.5)*dx) * (u[i+1,j] - u[i,j]) - 
					alpha((i-0.5)*dx) * (u[i,j] - u[i-1,j])) + 
					2*u[i,j] - u[i,j-1]
					, 10)

		print(i, j)

#		if u[i,j] != 0:
#			print(x[i], t[j], u[i,j])


fig = plt.figure()


#Writer = animation.writers['ffmpeg'] 
#writer = Writer(fps=35, metadata=dict(artist='Me'), bitrate=1800) 


ax = plt.axes(xlim = (3.45, 4.75), ylim=(-1.5, 1.5))
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
ax = plt.axes(projection='3d', xlim=(0, T), ylim=(2,5), zlim=(-0.5,1.5))
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


ax.set_title('SDC ex. 1 with time as space variable')
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('u')
ax.view_init(35, 180)


anim = animation.FuncAnimation(fig, animate, frames=len(t), 
		interval=1)
'''

#anim.save('1d_sdc_anim.mp4', writer=writer)

plt.show()





















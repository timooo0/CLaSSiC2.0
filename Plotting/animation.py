import helper
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

param, x, y, z = helper.getData()
fNum  = 0
atom = 0
distance = 1
structure = "triangle"

fig = plt.figure()
if structure == "line":
    ax = plt.axes(xlim=(-0.002, distance*(param[fNum]["atoms"])+0.002), ylim=(-distance*(param[fNum]["atoms"])/4-0.002, distance*(param[fNum]["atoms"])/4+0.002))
elif structure == "triangle":
     ax = plt.axes(xlim=(-1, 2*np.sqrt(param[fNum]["atoms"])+1), ylim=(-1, 2*np.sqrt(param[fNum]["atoms"])+1))
line, = ax.plot([], [], lw=2,color="b")
lines = []

for index in range(int(param[fNum]["atoms"])):
    lobj = ax.plot([],[],lw=2,color="b")[0]
    lines.append(lobj)

def init():
    for line in lines:
        line.set_data([], [])
    return lines

def animate(i):
    ax.set_title(f"time: {i} ps")
    for atom, line in enumerate(lines):
        if atom%1==0:
            if structure == "line":
                line.set_data([atom*distance, atom*distance+x[fNum, atom, i]], [0, y[fNum, atom, i]])
            elif structure == "triangle":
                sidelength = int(np.sqrt(param[fNum]["atoms"]))
                line.set_data([atom%sidelength*2 + atom/sidelength%2, atom%sidelength*2+ atom/sidelength%2+x[fNum, atom, i]], [atom/sidelength*2, atom/sidelength*2+y[fNum, atom, i]])
                # j+offset, 0.5*np.sqrt(3)*i, 0
            # if atom==0:
                # print(x[fNum, atom, i], y[fNum, atom, i])
    return lines

ani = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=range(0, x.shape[2], 1), interval=1000, blit=True)
# ani.save('matplot003.gif', writer='pillow')
plt.show()

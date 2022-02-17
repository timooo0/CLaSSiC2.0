import helper
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

fNum  = 0
param, x, y, z = helper.getData()
positions = helper.getPositions(fNum)
# print(positions.shape)
rescale = 0.5/np.max(x[fNum, 0, 0])

fig = plt.figure()

if param[fNum]["nDimensions"] == 1:
    ax = plt.axes(xlim=(-1, param[fNum]["nUnitCells"]), ylim=(-(param[fNum]["nUnitCells"])/2, (param[fNum]["nUnitCells"])/2))
    ax.axis("equal")
elif param[fNum]["nDimensions"] == 2:
    ax = plt.axes(xlim=(-1, param[fNum]["nUnitCells"]), ylim=(-1, param[fNum]["nUnitCells"]))
    ax.axis("equal")


lines = []
for index in range(int(param[fNum]["atoms"])):
    lobj = ax.plot([],[],"b-o", lw=2)[0]
    lines.append(lobj)

def init():
    for line in lines:
        line.set_data([], [])
    return lines

def animate(i):
    ax.set_title(f"time: {i} ps")
    for atom, line in enumerate(lines):
        line.set_data([positions[atom][0], positions[atom][0]+rescale*x[fNum, atom, i]], [positions[atom][1], positions[atom][1]+rescale*y[fNum, atom, i]])
    return lines

ani = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=range(0, x.shape[2], 1), interval=10, blit=True)
# ani.save('matplot003.gif', writer='pillow')
plt.show()

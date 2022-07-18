import numpy as np
import matplotlib.pyplot as plt

radius = 0.015
D = 2. * radius
y0 = 0.0
r_inlet_start = 4. * radius
r_inlet_end = 4. * radius
height = 10. * D

r = r_inlet_start
r2 = r_inlet_end
xlo = -height
xhi = 0.
xdiv = -4.0 * r

poly = [ (xlo, y0+r), (xdiv, y0+r), (xhi, y0+r2), (xhi, y0-r2), (xdiv, y0-r), (xlo, y0-r), (xlo, y0+r) ]

x = [el[0] for el in poly]
y = [el[1] for el in poly]

fig, ax = plt.subplots()

ax.plot(x,y, "x-" )

[ ax.annotate(str(i), xy=point, xycoords=ax.transData) for i, point in enumerate(poly[:-1]) ]

fig.savefig("temp.png")
fig.clf()
plt.close(fig)

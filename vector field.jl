using Plots

# create x and y arrays
x = range(-5, stop=5, length=10)
y = range(-5, stop=5, length=10)

# create arrays for the x and y components of the arrows
u = zeros(length(x))
v = -ones(length(y))

# create the quiver plot
quiver(x, y, (u, v))

savefig("vf.png")
### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 72906394-a47b-4e7a-98b4-e3c41e68dd4d
md"""# Tackpoint.jl
by *Luca Visser*
"""

# ╔═╡ 12fbea7a-e6bf-46a4-ab03-44bd10e7e256
md""" ## Introduction
This project focuses on developing algorithms that address the issue of optimizing a sailboat’s trajectory when a starting point and destination are given alongside static wind conditions. The underlying physics that govern the optimal path of a sailboat for a given set of conditions are highly coupled and dynamic, rendering the course very unintuitive to determine. Algorithms that are able to produce the path plan that takes the minimum amount of time to complete the course can be very helpful and the resulting work is presented in this report.

Algorithms developed in this project use the idea of calculating the Velocity Made Good (VMG) as a parameter relating the state of the sailboat at any given time to the time it would take to complete a given course.
Before going into details about how the models work, some basic sailing theory is introduced.

It is not possible for boats to sail directly into the wind, requiring the course of the boat to alternate between headings. This process is called "tacking" and is used commonly by sailors to make their way to a mark that is upwind. On a tack, the sailor will generally point the sailboat as close to the wind as possible while still keeping the winds blowing through the sails in a manner that provides aerodynamic lift to propel the boat.

Then the boat is turned away from the wind in slight increments in order to generate more forward lift on the sails allowing it to move with greater speed, but less directly toward the destination. This may be seen in the schematic diagram depicted in Figure 1.

The range of heading that does not produce any significant lift is called the no-go zone. It is indicated in red in Figure 2.
Figure 2: Points of sail – Red indicates no-go zone
"""


# ╔═╡ bbbf78d6-7bf0-4eb0-883b-9bc8f98a18d3
md""" ## Model
The VMG can be calculated by using the following expression where Vtrue is the velocity of the sailboat with respect to stationary ground and 𝜃 is the angle between current heading and the direction to destination. 𝑉𝑀𝐺=𝑉𝑡𝑟𝑢𝑒∗ 𝑐𝑜𝑠𝜃
However, the relationship between the true velocity of the sailboat and the true wind velocity depends on what assumptions are made. Some instances take into consideration the fact that a sailboat’s Speed Over Ground increases or decreases relative to the wind direction. In theory, a sailboat’s speed increases while sailing from upwind to a downwind direction. Other parameters to factor in include the sail boat’s specifics that depend on the make and design of the boat. These are the 'Velocity Increase Constant', ‘No – Go Zone’ and ‘Degree Interval’, which are normally provided by the manufacturer. The physical meaning of this constant expresses 'the sail boat increases speed by 10% for every X degrees from the wind'. Considering these factors, a True VMG can be calculated using the equations:
For Upwind:
For Downwind:
4
Hemanth Sarabu Optimizing Sailing Trajectories
Where,
𝑉𝑀𝐺=𝑉𝑒𝑙𝑜𝑐𝑖𝑡𝑦 𝑀𝑎𝑑𝑒 𝐺𝑜𝑜𝑑 𝑡𝑜𝑤𝑎𝑟𝑑𝑠 𝑑𝑒𝑠𝑡𝑖𝑛𝑎𝑡𝑖𝑜𝑛
𝑉𝑤=𝑉𝑒𝑙𝑜𝑐𝑖𝑡𝑦 𝑜𝑓 𝑊𝑖𝑛𝑑
𝜃𝑠=𝐴𝑛𝑔𝑙𝑒 𝑏𝑒𝑡𝑤𝑒𝑒𝑛 ℎ𝑒𝑎𝑑𝑖𝑛𝑔 𝑎𝑛𝑑 𝑑𝑒𝑠𝑡𝑖𝑛𝑎𝑡𝑖𝑜𝑛
𝜃0=𝑁𝑜−𝑔𝑜 𝑧𝑜𝑛𝑒
𝜃𝑦=𝐴𝑛𝑔𝑙𝑒 𝑏𝑒𝑡𝑤𝑒𝑒𝑛 𝑤𝑖𝑛𝑑 𝑑𝑖𝑟𝑒𝑐𝑡𝑖𝑜𝑛 𝑎𝑛𝑑 ℎ𝑒𝑎𝑑𝑖𝑛𝑔 β=𝑉𝑒𝑙𝑜𝑐𝑖𝑡𝑦 𝐼𝑛𝑐𝑟𝑒𝑎𝑠𝑒 𝐶𝑜𝑛𝑠𝑡𝑎𝑛𝑡 𝑖= 𝐷𝑒𝑔𝑟𝑒𝑒 𝐼𝑛𝑡𝑒𝑟𝑣𝑎𝑙

In the expressions shown above the no-go zone (𝜃𝑜) , Degree Interval (i) and Velocity Increase Constant (β) are usually provided by the manufacturer. As these are specific to the boat’s design, commonly quoted values are used throughout this project.
"""

# ╔═╡ 43c22300-f7ae-11ed-360f-c5cfdfef8406


# ╔═╡ 0952f6cf-cfe5-44df-9783-1ba7da0415bb


# ╔═╡ Cell order:
# ╟─72906394-a47b-4e7a-98b4-e3c41e68dd4d
# ╠═12fbea7a-e6bf-46a4-ab03-44bd10e7e256
# ╠═bbbf78d6-7bf0-4eb0-883b-9bc8f98a18d3
# ╠═43c22300-f7ae-11ed-360f-c5cfdfef8406
# ╠═0952f6cf-cfe5-44df-9783-1ba7da0415bb

### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 72906394-a47b-4e7a-98b4-e3c41e68dd4d
begin

using LinearAlgebra, Dates, Plots, PlutoUI, ForwardDiff

md"""# Tackpoint.jl
by *Luca Visser*
"""

end

# ╔═╡ 12fbea7a-e6bf-46a4-ab03-44bd10e7e256
begin
	md""" ## Introduction
	In the last decades, the sport of sailing has experienced an increasing impact
	of new technologies, and notably of scientific computing. Among all computational problems relevant for sailing, we are interested here in route planning
	and race strategy, i.e., the optimization of the yacht route.
	
	This project focuses on developing algorithms that address the issue of optimizing a sailboat’s trajectory when a starting point and destination are given alongside static wind conditions. The underlying physics that govern the optimal path of a sailboat for a given set of conditions are highly coupled and dynamic, rendering the course very unintuitive to determine. Algorithms that are able to produce the path plan that takes the minimum amount of time to complete the course can be very helpful.
	
	Algorithms developed in this project use the idea of calculating the Velocity Made Good (VMG) as a parameter relating the state of the sailboat at any given time to the time it would take to complete a given course. 
	
	Before going into details about how the models work, some basic sailing theory is introduced. It is not possible for boats to sail directly into the wind, requiring the course of the boat to alternate between headings. This process is called "tacking" and is used commonly by sailors to make their way to a mark that is upwind. On a tack, the sailor will generally point the sailboat as close to the wind as possible while still keeping the winds blowing through the sails in a manner that provides aerodynamic lift to propel the boat.
	
	Then the boat is turned away from the wind in slight increments in order to generate more forward lift on the sails, allowing it to move with greater speed, but less directly toward the destination. The range of heading that does not produce any significant lift is called the no-go zone.

	$(Resource("https://github.com/Lucmeister55/STMO-assignment/blob/main/images/sailing_intro.png?raw=true", :width => 1000, :align => "middle"))
	"""
end

# ╔═╡ 7c54bb74-2116-4540-9dbd-f10d49e3a3bd
md""" ## Methods
"""

# ╔═╡ bbbf78d6-7bf0-4eb0-883b-9bc8f98a18d3
md""" ### VMG
The VMG can be calculated by using the following expression where $V_{true}$ is the velocity of the sailboat with respect to stationary ground and $\theta_s$ is the angle between current heading and the direction to destination. 

``VMG=V_{true}∗\cos(\theta_s)``

However, the relationship between the true velocity of the sailboat and the true wind velocity depends on what assumptions are made. Some instances take into consideration the fact that a sailboat’s Speed Over Ground increases or decreases relative to the wind direction. In theory, a sailboat’s speed increases while sailing from upwind to a downwind direction. Other parameters to factor in include the sail boat’s specifics that depend on the make and design of the boat. These are the 'Velocity Increase Constant', ‘No – Go Zone’ and ‘Degree Interval’, which are normally provided by the manufacturer. The physical meaning of this constant expresses 'the sail boat increases speed by 10% for every X degrees from the wind'. Considering these factors, a True VMG can be calculated using the equations:

For Upwind: 

$$VMG=\frac{V_w}{\cos(\theta_0)}*(1+\beta)^{\frac{|\theta_0-\theta_\gamma|}{i}}*\cos(\theta_s)$$

For Downwind: 

$$VMG=\frac{V_w}{\cos(\theta_0)}*(1+\beta)^{\frac{|180°-\theta_0-\theta_\gamma|}{i}}*\cos(\theta_s)$$

VMG = Velocity Made Good towards destination\
$V_w$ = Velocity of wind\
$\theta_s$ = Angle between heading and destination\
$\theta_0$ = No-go zone\
$\theta_\gamma$ = Angle between wind direction and heading\
$\beta$ = Velocity increase constant\
i = Degree Interval

In the expressions shown above the no-go zone ($\theta_0$), Degree Interval (i) and Velocity Increase Constant ($\beta$) are usually provided by the manufacturer. As these are specific to the boat’s design, commonly quoted values are used throughout this project.
"""

# ╔═╡ cea4baae-51e5-43a6-a64f-a11e3d474fb6
md""" ### Assumptions
In order to simplify the optimization problem, the following assumptions are made and incorporated into the models:
- True velocity of boat at any point of time is a function of the velocity of wind and heading relative to wind direction only. This implies that at any point these are the only two variables required to calculate how fast the boat will be moving with respect to ground (true velocity) in the direction it is moving.
- Effects of momentum are not considered. As the expressions shown earlier allow for the calculation of velocity instantaneously, the time-dynamic effects are not modelled. This means that acceleration and deceleration are not taken into account and if the boat moves from a current position to the next position, the momentum is not carried or lost but the position’s assigned momentum is taken up.
- Effects of drag are not considered. The effects of drag considered are only those specified by the velocity increase constant and the no-go zone, both provided by the manufacturer. The effects of drag due to the wind and the water due to the form of the boat are not taken into consideration. In real life however, this phenomenon will have significant effects limiting the maximum velocity of the sailboat.
- Effects of manoeuvres on the momentum of the sailboat are not considered. For example, tacking causes the sailboat to lose momentum due to the associated drag and loss of lift during the procedure. The loss of momentum is a transient process and is not modelled for this project.
"""

# ╔═╡ 1f442a09-9b34-4644-949d-b649187470c2
md"""### Decision Variables
The objective of the algorithms is to produce paths that take the least amount of time to complete a given course. In order to describe the path taken by the sailboats, the variables required are the x and y coordinates of the sailboat at a given time. Instead of determining the coordinates for each time step, 
they are determined in consecutive steps-this means the coordinates $X^{k+1}$ are determined relative to $X^k$ in distance and direction as opposed to determining $X^{k+1}$ as a function of time step $t^{k+1}$.

The problem of producing a locus of X, Y coordinates that are successive in nature can be transformed into producing a direction vector and a distance vector from each of the points to the next point. An example of this is- if $X^k$ is known, in order to determine the next position $X^{k+1}$ a direction vector and a step size are sufficient and required with relation to $X^k$. Hence, the modelling problem is reduced to finding the optimal sequence of direction vectors (heading) if the step size is set to be constant. 

In other standard optimization problems tackled during the course, the design space did not include a time dimension. However, this problem indeed has a time dimension and the variable (heading) will have to be determined as the time marches forward (in distance steps and not time steps).
"""

# ╔═╡ 75802d02-ee43-4028-ab28-483f9af7a54c
md"""### Objective
The objective of the algorithm is to find the path of least resistance, or the path that takes the least amount of time to complete a given course. The time taken is calculated in a discretized manner: The velocity of the sailboat is found for each step and is used to divide the step size.

$\textrm{time step k} = \frac{\textrm{step size}}{\textrm{velocity in direction of step}}$

The time steps are added for the duration of the course in order to register the total time taken.

$\textrm{time elapsed} = \sum{\textrm{time step k}}$

The objective is to optimize the path to achieve the lowest time elapsed value for any given course.
"""

# ╔═╡ c279c555-7d66-481e-b332-8519fba0f950
md"""### Constraints
This path optimization problem only has one constraint – at no instance should the boat bear a heading into the no-go zone. In real life however, during a tack, the boat indeed faces the no-go zone briefly before regaining lift in the sails. It is its momentum that allows the boat to steer away from the no-go zone and into a heading that permits forward motion. In the models presented in the report, it is 
important to remember that the boat’s momentum terms are not incorporated. This means that if at any point the algorithm produces a path that involves a heading into the no-go zone, it loses its velocity and cannot make further steps. The algorithms get past this flaw by limiting the heading from entering the no-go zone at all times.

When travelling upwind,

$\theta_w > \theta_0$\

$\theta_0$ = no-go zone\
$\theta_w$ = heading relative to wind direction\
"""

# ╔═╡ 3fdab425-52ef-4be5-a8a2-425cd2a786ef
md"""### Model
The Single Tack Method (STM) is developed so that the path produced by the model is extremely easy to navigate. This is achieved by implicitly allowing the algorithm to produce a path between two points that incorporates a maximum of one tack (i.e. two leg journey). The algorithm is then required to produce the tack length and angle that results in the minimum time elapsed value. It is important to note that the manner in which ease of navigation is incorporated into this model is by limiting the number of tacks. By increasing the number of tacks, the solution process can get very complex. 

To explain this, a scenario is discussed – If arbitrary location A is considered the starting point and a destination is set at B, to produce a path that can only have a single tack, the problem is reduced to finding a point Ts (tack point) such that when the boat travels from A to Ts, Ts to B, it takes less amount 
of time than if the boat sailed from A to T, T to B. Here T is an arbitrary tack location. If two tacks were allowed, this would mean that the algorithm would have to find the best combination of two points that result in the smallest time elapsed value. Due to the highly coupled nature of the problem, navigating through two dimensions to find the best combination would be very difficult, and hence is out of the scope of this project.
"""

# ╔═╡ d0019978-aa34-4f54-b4b7-4d607274b5a6
md"""## Implementation
"""

# ╔═╡ 913610b3-1ecd-49e1-a692-41c86ef0431e
md""" ### Objective function
"""

# ╔═╡ f65b308e-8764-4c2b-851b-aa16d1277041
function my_acos(x)
    if isapprox(x, 1)
        return zero(x)
    elseif isapprox(x, -1)
        return one(x) * π
    else
        return acos(x)
    end
end

# ╔═╡ 55281d8b-38a7-465c-8419-e90ac5c316ab
function pathtime(xt, cons, theta_nogo = deg2rad(40))
	x0, x_tar_c, vel_wind, wind_dir, maxiter = cons
	
    # x_tar_c = active target
    # x0 = current position
    # xt = tack coordinates

    velcons = 3 # velocity increase constant
    deg_int = 5 # constant provided by manufacture
	step_inc = 0.1
    
    # tack vector 1
    tack_vect_1 = xt - x0 
    # tack vector 2
    tack_vect_2 = x_tar_c - xt
    
    # divide it into increments of .1 metres
    steps_1 = norm(xt - x0)/step_inc
    # set initial step
    step_c_1 = 0
    heading = tack_vect_1
    xc = x0
    time_elapsed_1 = 0
    
    while step_c_1 < steps_1
        # vector to target from current x
        target_vect = x_tar_c - xc 
        
        # calculate angle between heading and target vector
        theta_rd = real(my_acos(dot(heading, target_vect)/(norm(heading)*norm(target_vect))))
        
        # calculate angle between heading and wind direction
        theta_rw = real(my_acos(dot(heading,wind_dir)/(norm(heading)*norm(wind_dir))))
        
        if dot(wind_dir, heading)/(norm(wind_dir)*norm(heading)) < 0
            # calculate velocity made good
            u = cos(theta_rd)*(vel_wind/cos(theta_nogo))*(1+0.01*velcons)^((abs(theta_rw-theta_nogo)*180)/(pi*deg_int))
            # check if within deadzone
            v = (my_acos(dot(heading,-wind_dir)/(norm(heading)*norm(-wind_dir)))>theta_nogo)*u
        elseif dot(wind_dir,heading)/(norm(wind_dir)*norm(heading)) > 0
            # calculate velocity made good if downwind
            v = cos(theta_rd)*(vel_wind/cos(theta_nogo))*(1+0.01*velcons)^((abs(pi-theta_rw-theta_nogo)*180)/(pi*deg_int))
        end
        
        if abs(cos(theta_rd)) > 1e-6
		    vel_tack = v / cos(theta_rd)
		else
		    vel_tack = 0.0  # or any other appropriate value
		end
        
        if vel_tack > 10e-2
            # time elapsed in increment
            time_increment_1 = 0.1/vel_tack 
            time_elapsed_1 += time_increment_1
            xc += 0.1*heading/norm(heading)
            step_c_1 += 1
        elseif vel_tack < 10e-2
            time_elapsed_1 = 10e6
            break
        end
    end
    
    t_1 = time_elapsed_1
    
    # divide it into increments of .1 metres
    steps_2 = norm(x_tar_c - xt)/step_inc
    # set initial step
    step_c_2 = 0
    heading = tack_vect_2
    xc = xt
    time_elapsed_2 = 0
    while step_c_2 < steps_2
        target_vect = x_tar_c - xc
        theta_rd = real(my_acos(dot(heading,target_vect)/(norm(heading)*norm(target_vect))))
        theta_rw = real(my_acos(dot(heading,wind_dir)/(norm(heading)*norm(wind_dir))))
        if dot(wind_dir,heading)/(norm(wind_dir)*norm(heading)) < 0
            u = cos(theta_rd)*(vel_wind/cos(theta_nogo))*(1+0.01*velcons)^((abs(theta_rw-theta_nogo)*180)/(pi*deg_int))
            v = (my_acos(dot(heading,-wind_dir)/(norm(heading)*norm(-wind_dir)))>theta_nogo)*u
        elseif dot(wind_dir,heading)/(norm(wind_dir)*norm(heading)) > 0
            v = cos(theta_rd)*(vel_wind/cos(theta_nogo))*(1+0.01*velcons)^((abs(pi-theta_rw-theta_nogo)*180)/(pi*deg_int))
        end
		
		if abs(cos(theta_rd)) > 1e-6
		    vel_tack = v / cos(theta_rd)
		else
		    vel_tack = 0.0  # or any other appropriate value
		end
		
        if vel_tack > 10e-2
            time_increment_2 = 0.1/vel_tack
            time_elapsed_2 = time_elapsed_2 + time_increment_2
            xc = 0.1*heading/norm(heading) + xc
            step_c_2 = step_c_2 + 1
        elseif vel_tack < 10e-2
            time_elapsed_2 = 10e6
            break
        end
    end
    t_2 = time_elapsed_2
    t = t_1+t_2
    return t
end

# ╔═╡ 76c07283-c73b-4b3d-9396-e7f704bf7312
function grad_f(func, point, constants)
    f(x) = func(x, constants)
    gradient = ForwardDiff.gradient(f, point)
    return gradient
end

# ╔═╡ 99f4fcca-19bc-4a80-8472-198ed377abd0
function hess_f(func, point, constants)
    f(x) = func(x, constants)
    hessian = ForwardDiff.hessian(f, point)
    return hessian
end

# ╔═╡ 81c03358-1a48-4c33-802a-f2bfbf89f6d7
md""" ### Adaptive Step Size
"""

# ╔═╡ 0cbef7ea-c480-4679-88a2-f9e929f3c676
md""" #### Backtracking Line Search
"""

# ╔═╡ b1b865d3-e4cd-4850-adba-9462b5930929
function backtracking_line_search(f, x0, Dx, cons, c=0.01, rho=0.1)
    
	alpha=30
	
    while (f(x0 + alpha*Dx, cons) > f(x0, cons) + alpha * c * transpose(-Dx)*Dx)
        alpha *= rho
	end
	
    return alpha
end

# ╔═╡ 0d55edfd-f003-466a-a12c-b4463f28168b
md""" ### Constraints
"""

# ╔═╡ bfb209de-c20b-455a-b45e-9f544d4877b0
function constraint_satisfied(xt, cons, theta_nogo = deg2rad(50))
    x0, x_tar_c, vel_wind, wind_dir, maxiter = cons
    
    # Calculate the angle between wind direction and tack vectors
    tack_vect_1 = xt - x0
    tack_vect_2 = x_tar_c - xt
    
    theta_rw_1 = real(my_acos(dot(tack_vect_1, -wind_dir) / (norm(tack_vect_1) * norm(-wind_dir))))
    theta_rw_2 = real(my_acos(dot(tack_vect_2, -wind_dir) / (norm(tack_vect_2) * norm(-wind_dir))))
    
    # Check if both angle constraints are satisfied
    if theta_rw_1 >= theta_nogo && theta_rw_2 >= theta_nogo
        return true
    else
        return false
    end
end

# ╔═╡ 159b8d41-a5ce-494a-9566-df34547b1a8d
md""" ### Tackpoint Estimation
"""

# ╔═╡ 7700b89e-8322-4d42-9b35-dd391731b5e3
md""" ### Initial Estimate
"""

# ╔═╡ 17dd25fd-1640-4f4b-9349-e35d73a897b8
function initial_xt_dev(x0, x_tar_c)
	#calculating initial xt---------------------------------------------------------
    vect_tar = (x_tar_c - x0) #generate vector to active target
    mid = 0.5 .* vect_tar + x0 #place point in the middle
    deviation = ((x_tar_c[2]-x0[2])^2/(x_tar_c[1]-x0[1]))/(atan(deg2rad(15))) #define deviation of point from vect_tar
    xt0 = [-deviation, 0] + mid #initial tack point
	return xt0
end

# ╔═╡ 94e8d30c-9280-427c-a5bf-c3370712e69f
function calculate_y(x, point, vector)
    x_1, y_1 = point
    a, b = vector
    
    # Calculate the corresponding y value using the equation of the line
    y = (a * x_1 + b * y_1 - a * x) / b
    
    return y
end

# ╔═╡ d8f88426-d2a1-4ae1-8571-2ea80975a4fd
function initial_xt_rand(x0, x_tar_c, wind_dir, theta_nogo = deg2rad(50))

	vect_tar = (x_tar_c - x0) #generate vector to active target
    mid = 0.5 .* vect_tar + x0 #place point in the middle
	
    max_iter = 1000  # Maximum number of attempts to find a feasible starting point
    
    for i in 1:max_iter
        xt = [rand(-100:100), rand(round(calculate_y(-100, mid, wind_dir)):round(calculate_y(100, mid, wind_dir)))]  # Generate a random point
        
        # Check if the angle constraint is satisfied
        tack_vect_1 = xt - x0
		tack_vect_2 = x_tar_c - xt
		
        theta_rw_1 = my_acos(dot(tack_vect_1, -wind_dir) / (norm(tack_vect_1) * norm(-wind_dir)))
		theta_rw_2 = my_acos(dot(tack_vect_2, -wind_dir) / (norm(tack_vect_2) * norm(-wind_dir)))
        
        if theta_rw_1 >= theta_nogo && theta_rw_2 >= theta_nogo
            return theta_rw_1, theta_rw_2, xt
        end
    end

	error("Unable to find a feasible starting point within the given number of iterations.")
end


# ╔═╡ e9c2610f-f94d-4426-be0f-d11d8e1498ab
md""" #### Line Search (Brute Force)
The starting point and the current destination is specified in this routine. The programme then calculates an initial feasible point for the tack location. The model employs a pattern search algorithm incorporating an accelerating/decelerating step size. The pattern search algorithm requires to start in the feasible region because starting in the no-go zone results in an infinite path time value.

Once, the pattern search algorithm is initialised, the path time is calculated in each of the probe directions (x+,y+,x-,y-) using the function handle pathtime. This function handle accepts the starting point of the sailboat, the tack point and the destination point and integrates along the two legs of the journey to find the time elapsed (or the path time). This is the objective function value that needs to be minimised. 

With the objective function values from all four probes, the pattern search algorithm chooses the directions in x and y that favour the minimizing of the path time and determines the new tack point. Depending on the search directions (in X and Y) recorded, the pattern search algorithm updates the acceleration terms to reduce the number of pattern moves required to achieve convergence. Once convergence is reached, the optimal tack point (xt) and the optimal time elapsed values are returned.
"""

# ╔═╡ 546844bf-3e36-4e9f-9bc9-aeb92649b108
function tackpoint_LS(f, cons, xt0)
	x0, x_tar_c, vel_wind, wind_dir, maxiter = cons
	
	xt = xt0

    iter = 0

    #Introducing probes for pattern search
    plength_x = [10,0] #probe length in x
    plength_y = [0,10] #probe length in y
    ax = 0.1 #initial stepsize
    ay = 0.1 #initial stepsize
    adpxt = 1 #initialise acceleration determining parameter in x
    adpyt = 1 #initialise accleration determining parameter in y

    while iter < maxiter # set maximum number of iterations
		
        # calculate path times at probes
        pathtime_xtc = f(xt, cons) # path time if current tack point is used
        pathtime_xtc_xplus = f(xt + plength_x, cons)
        pathtime_xtc_xminus = f(xt - plength_x, cons)
        pathtime_xtc_yplus = f(xt + plength_y, cons)
        pathtime_xtc_yminus = f(xt - plength_y, cons)

        # probing in x and updating xt
        # find minimum path time
        if min(pathtime_xtc, pathtime_xtc_xplus, pathtime_xtc_xminus) == pathtime_xtc
            xt = xt
            sx = 0 # set search direction
            adpx = 0 # set acceleration term
        elseif min(pathtime_xtc, pathtime_xtc_xplus, pathtime_xtc_xminus) == pathtime_xtc_xplus
            sx = 1 # set search direction
            adpx = 1 # set acceleration term
        elseif min(pathtime_xtc, pathtime_xtc_xplus, pathtime_xtc_xminus) == pathtime_xtc_xminus
            sx = -1 # set search direction
            adpx = -1 # set acceleration term
        end

        # probing in y and updating xt
        if min(pathtime_xtc, pathtime_xtc_yplus, pathtime_xtc_yminus) == pathtime_xtc
            sy = 0 # set search direction
            adpy = 0 # set acceleration term
        elseif min(pathtime_xtc, pathtime_xtc_yplus, pathtime_xtc_yminus) == pathtime_xtc_yplus
            sy = 1 # set search direction
            adpy = 1 # set acceleration term
        elseif min(pathtime_xtc, pathtime_xtc_yplus, pathtime_xtc_yminus) == pathtime_xtc_yminus
            sy = -1 # set search direction
            adpy = -1 # set acceleration term
        end
        iter += 1
        ax = ax*1.15^(adpx*adpxt) # update acceleration in x
        ay = ay*1.15^(adpy*adpyt) # update acceleration in y
        xt_p = xt # update xt previous
        xt = xt + [ax*sx, ay*sy] # update tack position
        adpxt = adpx # store acceleration terms
        adpyt = adpy
        minpathtime = f(xt, cons) # calculate corresponding minpathtime
        minpathtime_p = f(xt_p, cons)

        # prevent entering into infeasible region
        if minpathtime_p < minpathtime
            xt = xt_p
            minpathtime = minpathtime_p
            ax = 0.1 # reset acceleration
            ay = 0.1 # reset acceleration
        end

		println("Iteration: " * string(iter))
		println("New tack point: " * string(xt))
		println()
    end
    return xt
end

# ╔═╡ d628a011-a1fb-4947-8da1-2a1c7e2cf048
md""" #### Gradient Descent
"""

# ╔═╡ 4622c7a9-a35c-44ef-9c4a-e6d09faa6e08
function tackpoint_GD(f, cons, xt0, nu=1e-3)
	
	x0, x_tar_c, vel_wind, wind_dir, maxiter = cons
	
	xt = xt0

    iter = 0

    while iter < maxiter # set maximum number of iterations
        
		Dx = -grad_f(f, xt, cons)

		if norm(-Dx) <= nu
            break  # converged
		end
		
		t = backtracking_line_search(f, xt, Dx, cons)

		if constraint_satisfied(xt+Dx*t, cons)
            xt += t * Dx
        else
            xt -= t * Dx
        end

        iter += 1

		println("Iteration: " * string(iter))
		println("Search direction: " * string(Dx))
		println("Step size: " * string(t))
        println("New tack point: " * string(xt))
		println()
    end
    return xt
end

# ╔═╡ 578f8606-4873-4c96-8112-1199fa486af1
md""" #### Newton's Method
"""

# ╔═╡ f0a5869c-e3eb-425d-aab4-02b8018ee6a3
function tackpoint_NM(f, cons, xt0 = nothing, epsilon=1e-3)
	
	x0, x_tar_c, vel_wind, wind_dir, maxiter = cons
	
	xt = xt0

    iter = 0

    while iter < maxiter # set maximum number of iterations
        
		Dx = -hess_f(f, xt, cons) \ grad_f(f, xt, cons)

		println(hess_f(f, xt, cons))
		println(grad_f(f, xt, cons))
		println(Dx)
		println(-dot(grad_f(f, xt, cons), Dx) / 2)
		
        if -dot(grad_f(f, xt, cons), Dx) / 2 <= epsilon  # stopping criterion
            break  # converged
        end
		
		t = backtracking_line_search(f, xt, Dx, cons)

		if constraint(xt, cons)
            xt += t * Dx
        else
            xt -= t * Dx
        end

        iter += 1

		println("Iteration: " * string(iter))
		println("Search direction: " * string(Dx))
		println("Step size: " * string(t))
        println("New tack point: " * string(xt))
		println()
    end
    return xt
end

# ╔═╡ 6ff92112-41d4-47bb-9eb0-4a2499f386fd
md""" ### Calculate Path
Function to calculate the path produced by this algorithm.
"""

# ╔═╡ c315d2b8-1a38-4ae1-b6d1-151b1a99cd20
function track(x0, x_tar, wind_dir, target_tol)
    xp = x0 # current position
    # Initialize target tracking
    num_targets = size(x_tar, 1) # get number of targets
    target_tracker = zeros(Int64, num_targets) # create file indicating no targets have been achieved
    itt = 1 # target tracker iteration
    while itt < num_targets + 1
        target_tracker[itt] = 0
        itt += 1
    end

	next = "tackpoint"

    x_tar_c = [x_tar[1,1], x_tar[1,2]] # set initial target
    tar_act = 1 # active target is 1

	# Compute and print headings relative to wind
	v_tw = x_tar_c - xp 
	a_tw = real(my_acos(dot(-wind_dir, v_tw)/(norm(-wind_dir)*norm(v_tw))))

	tacknr = 1
	
	println()
	println("Heading relative to wind before tack "*string(tacknr)*": " *string(Int64(round(rad2deg(a_tw))))*"°")

    # Target tracking loop
    x = zeros(0)
    y = zeros(0)
	
    while norm(x_tar_c - xp) > target_tol / 10 # if target is not achieved
        step_direction = (x_tar_c - xp) / norm(x_tar_c - xp)
        stepsize = 1
        xp_p = xp # store previous location
        xp = step_direction .* stepsize .+ xp
        append!(x, xp[1])
        append!(y, xp[2])
        if norm(x_tar_c - xp) < target_tol
            if tar_act < size(x_tar, 1)
                target_tracker[tar_act] = 1
                x_tar_c = [x_tar[tar_act + 1, 1], x_tar[tar_act + 1, 2]]
                tar_act += 1
				next = "destination"
				
				v_tw = x_tar_c - xp 
				a_tw = real(my_acos(dot(-wind_dir, v_tw)/(norm(-wind_dir)*norm(v_tw))))
				tacknr += 1
				println("Heading relative to wind before tack "*string(tacknr)*": " *string(Int64(round(rad2deg(a_tw))))*"°")
            end
        end
    end
    return x, y
end

# ╔═╡ d4da9d8c-2064-4d1f-b397-4092ef0de660
md""" ### Plot Path
"""

# ╔═╡ 52a9641b-ea4c-45ce-b52f-ed3d47570a42
function plot_path(xt, cons, target_tol = 5)
	x0, x_tar_all, vel_wind, wind_dir, maxiter = cons

	x_tar = Array{Float64}(undef, 0, 2)

	for i in 1:length(x_tar_all)
		x_tar = [x_tar; xt[i][1] xt[i][2]; x_tar_all[i][1] x_tar_all[i][2]]
	end

	x, y = track(x0, x_tar, wind_dir, target_tol)

	# Plot route
	scatter(x, y, label = "Boat route")
	
	# Plot the wind vector
	xw = [x_tar[2, 1], x_tar[2, 1]+wind_dir[1]*10]
	yw = [x_tar[2, 2], x_tar[2, 2]+wind_dir[2]*10]
	plot!(xw, yw, arrow=true, arrowsize=0.5, label = "Wind direction")

	colors = ["green", "yellow", "orange", "red"]
	
	for i in 1:length(x_tar_all)
		scatter!([x_tar_all[i][1]], [x_tar_all[i][2]], color = colors[i], label = "Marker "*string(i), markersize = 10)
	end

	plot!()
end

# ╔═╡ 3fc06973-839c-47d4-b0ca-2cdef798d5ca
md""" ## Results
"""

# ╔═╡ 7c4e4208-944b-4c7e-a2fd-6275d2347a19
md""" ### Basic course
"""

# ╔═╡ fe1536be-8589-48ea-bbb8-9df739c4766c
md""" #### Parameter Options:
**Starting point**\
`X` : $(@bind x0_x Slider(-50:50, default=0, show_value=true))\
`Y` : $(@bind x0_y Slider(-50:50, default=0, show_value=true))\
**Destination point**\
`X` : $(@bind x_tar_x Slider(-50:50, default=20, show_value=true))\
`Y` : $(@bind x_tar_y Slider(-50:50, default=50, show_value=true))\
**Wind**\
`Wind velocity (m/s)` : $(@bind vel_wind Slider(1:10, default=5, show_value=true))\
`Wind direction (origin)` : $(@bind wind_dir_temp Select(["N", "NE", "E", "SE", "S", "SW", "W", "NW"], default="N"))\
**Algorithm**\
`Maximum iterations` : $(@bind maxiter Slider(100:100:1000, default=100, show_value=true))\
"""

# ╔═╡ 868a916b-1aa3-4f6d-acf8-84f9b955b3e5
begin
	if wind_dir_temp == "N"
		wind_dir = [0, -1]
	elseif wind_dir_temp == "NE"
		wind_dir = [-1, -1]
	elseif wind_dir_temp == "E"
		wind_dir = [-1, 0]
	elseif wind_dir_temp == "SE"
		wind_dir = [-1, 1]
	elseif wind_dir_temp == "S"
		wind_dir = [0, 1]
	elseif wind_dir_temp == "SW"
		wind_dir = [1, 1]
	elseif wind_dir_temp == "W"
		wind_dir = [1, 0]
	else
		wind_dir = [1, -1]
	end
	
	#Routine to determine the best tack point between two marks and calculate time
	x_tar_c = [x_tar_x, x_tar_y] #active destination
	x0 = [x0_x, x0_y] #starting point
end

# ╔═╡ 5c858f72-fcda-45b5-b9a3-036e5cc2d9e9
md""" #### Line Search
"""

# ╔═╡ 1885527b-2fb8-46e1-99a5-a21e0efd0358
let
	xt0 = initial_xt_dev(x0, x_tar_c)
	cons = [x0, x_tar_c, vel_wind, wind_dir, maxiter]
	
	t0 = now() # record start time of function
	xt = tackpoint_LS(pathtime, cons, xt0)

	# Record elapsed time
	solvertime = (now() - t0)
	println("Time until convergence: " * string(solvertime))
	println("Total pathtime: " * string(pathtime(xt, cons)))

	# Plot path
	cons = [x0, [x_tar_c], vel_wind, wind_dir, maxiter]
	plot_path([xt], cons)
end

# ╔═╡ cb3e7da3-4146-477e-8f24-6a55a1e800c7
md""" #### Gradient Descent
"""

# ╔═╡ fdc4ba68-07b7-4dcf-8e6b-d4e55a4f8826
begin
	theta_rw_1_GD, theta_rw_2_GD, xt0_GD = initial_xt_rand(x0, x_tar_c, wind_dir)
		println("Initial point: "*string(xt0_GD))
		println("Initial heading 1: "*string(rad2deg(theta_rw_1_GD)))
		println("Initial heading 2: "*string(rad2deg(theta_rw_2_GD)))
		println()
end

# ╔═╡ 664f2131-553d-4a89-ac5f-cf8138c59b4d
let
	cons = [x0, x_tar_c, vel_wind, wind_dir, maxiter]
	t0 = now() # record start time of function
	xt_GD = tackpoint_GD(pathtime, cons, xt0_GD)

	# Record elapsed time
	solvertime = (now() - t0)
	println("Time until convergence: " * string(solvertime))
	println("Total pathtime: " * string(pathtime(xt_GD, cons)))

	# Plot path
	cons_plot = [x0, [x_tar_c], vel_wind, wind_dir, maxiter]
	plot_path([xt_GD], cons_plot)
end

# ╔═╡ c8bb9329-3c83-4cff-83bc-de79ab29c5a4
md""" #### Newton's method
"""

# ╔═╡ f64da374-d99e-4be1-87ca-23261c56a26b
begin
	theta_rw_1_NM, theta_rw_2_NM, xt0_NM = initial_xt_rand(x0, x_tar_c, wind_dir)
		println("Initial point: "*string(xt0_NM))
		println("Initial heading 1: "*string(rad2deg(theta_rw_1_NM)))
		println("Initial heading 2: "*string(rad2deg(theta_rw_2_NM)))
		println()
end

# ╔═╡ 30cbcfc4-0780-47d2-aec0-d990436e1693
let
	cons = [x0, x_tar_c, vel_wind, wind_dir, maxiter]
	t0 = now() # record start time of function
	xt_NM = tackpoint_NM(pathtime, cons, xt0_NM)

	# Record elapsed time
	solvertime = (now() - t0)
	println("Time until convergence: " * string(solvertime))
	println("Total pathtime: " * string(pathtime(xt_NM, cons)))

	# Plot path
	cons_plot = [x0, [x_tar_c], vel_wind, wind_dir, maxiter]
	plot_path([xt_NM], cons_plot)
end

# ╔═╡ a1f98368-2ecb-40bc-80d2-b55d4224f2fc
md""" ### Racecourse
"""

# ╔═╡ e0b1b009-e7a4-47dd-bd17-e32377b2d566
md""" #### Parameter Options:
**Start/end point**\
`X` : $(@bind x0_x_rc Slider(-50:50, default=0, show_value=true))\
`Y` : $(@bind x0_y_rc Slider(-50:50, default=0, show_value=true))\
**Marker 1**\
`X` : $(@bind x_tar_x_rc1 Slider(-50:50, default=20, show_value=true))\
`Y` : $(@bind x_tar_y_rc1 Slider(-50:50, default=50, show_value=true))\
**Marker 2**\
`X` : $(@bind x_tar_x_rc2 Slider(-50:50, default=-50, show_value=true))\
`Y` : $(@bind x_tar_y_rc2 Slider(-50:50, default=0, show_value=true))\
**Marker 3**\
`X` : $(@bind x_tar_x_rc3 Slider(-50:50, default=20, show_value=true))\
`Y` : $(@bind x_tar_y_rc3 Slider(-50:50, default=-50, show_value=true))\
**Wind**\
`Wind velocity (m/s)` : $(@bind vel_wind_rc Slider(1:10, default=5, show_value=true))\
`Wind direction (origin)` : $(@bind wind_dir_temp_rc Select(["N", "NE", "E", "SE", "S", "SW", "W", "NW"], default="N"))\
**Algorithm**\
`Maximum iterations` : $(@bind maxiter_rc Slider(100:100:1000, default=100, show_value=true))\
"""

# ╔═╡ 13bdfa35-562b-49a5-978a-26c4a4920359
begin
	if wind_dir_temp_rc == "N"
		wind_dir_rc = [0, -1]
	elseif wind_dir_temp_rc == "NE"
		wind_dir_rc = [-1, -1]
	elseif wind_dir_temp_rc == "E"
		wind_dir_rc = [-1, 0]
	elseif wind_dir_temp_rc == "SE"
		wind_dir_rc = [-1, 1]
	elseif wind_dir_temp_rc == "S"
		wind_dir_rc = [0, 1]
	elseif wind_dir_temp_rc == "SW"
		wind_dir_rc = [1, 1]
	elseif wind_dir_temp_rc == "W"
		wind_dir_rc = [1, 0]
	else
		wind_dir_rc = [1, -1]
	end

	x0_rc = [x0_x_rc, x0_y_rc] #starting point
	x_tar_all = [x0_rc, [x_tar_x_rc1, x_tar_y_rc1], [x_tar_x_rc2, x_tar_y_rc2], [x_tar_x_rc3, x_tar_y_rc3], x0_rc]
end

# ╔═╡ f76d1f42-d5b9-40b5-9f31-667c04aee909
let
	xt_GD_rc = Vector{Float64}[]
	pt = 0
	for i in 1:length(x_tar_all)-1
		x0_rc = x_tar_all[i]
		x_tar_c_rc = x_tar_all[i+1]
		theta_rw_1_rc, theta_rw_2_rc, xt0_rc = initial_xt_rand(x0_rc, x_tar_c_rc, wind_dir_rc)
		println("Initial point: "*string(xt0_rc))
		println("Initial heading 1: "*string(rad2deg(theta_rw_1_rc)))
		println("Initial heading 2: "*string(rad2deg(theta_rw_2_rc)))
		println()
		cons_rc = [x0_rc, x_tar_c, vel_wind_rc, wind_dir_rc, maxiter_rc]
		
		t0 = now() # record start time of function
		xt_GD_temp = tackpoint_GD(pathtime, cons_rc, xt0_rc)
		push!(xt_GD_rc, xt_GD_temp)
		pt_temp = pathtime(xt_GD_temp, cons_rc)
		pt += pt_temp
	
		# Record elapsed time
		solvertime = (now() - t0)
		println()
		println("Time until convergence: " * string(solvertime))
		println("Total pathtime: " * string(pt_temp))
		println()
	end

	println()
	println("Racecourse pathtime: " * string(pt))

	# Plot path
	cons_rc = [x_tar_all[1], x_tar_all[2:end], vel_wind_rc, wind_dir_rc, maxiter_rc]
	plot_path(xt_GD_rc, cons_rc)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
ForwardDiff = "~0.10.35"
Plots = "~1.38.9"
PlutoUI = "~0.7.50"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "8e2970c58d4d38916e7fcf4bbbfbc3ad4fc7fcea"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "be6ab11021cd29f0344d5c4357b163af05a48cba"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.21.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "96d823b94ba8d187a6d8f0826e731195a74b90e9"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.0"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "d014972cd6f5afb1f8cd7adf000b7a966d62c304"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.5"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f670f269909a9114df1380cc0fcaa316fff655fb"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.5+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ba9eca9f8bdb787c6f3cf52cb4a404c0e349a0d1"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.5"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "099e356f267354f46ba65087981a77da23a279b7"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.0"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "a5aef8d4a6e8d81f171b2bd4be5265b01384c74c"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.10"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "d03ef538114b38f89d66776f2d8fdc0280f90621"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.12"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "259e206946c293698122f63e2b513a7c99a244e8"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─72906394-a47b-4e7a-98b4-e3c41e68dd4d
# ╟─12fbea7a-e6bf-46a4-ab03-44bd10e7e256
# ╟─7c54bb74-2116-4540-9dbd-f10d49e3a3bd
# ╟─bbbf78d6-7bf0-4eb0-883b-9bc8f98a18d3
# ╟─cea4baae-51e5-43a6-a64f-a11e3d474fb6
# ╟─1f442a09-9b34-4644-949d-b649187470c2
# ╟─75802d02-ee43-4028-ab28-483f9af7a54c
# ╟─c279c555-7d66-481e-b332-8519fba0f950
# ╟─3fdab425-52ef-4be5-a8a2-425cd2a786ef
# ╟─d0019978-aa34-4f54-b4b7-4d607274b5a6
# ╟─913610b3-1ecd-49e1-a692-41c86ef0431e
# ╠═f65b308e-8764-4c2b-851b-aa16d1277041
# ╠═55281d8b-38a7-465c-8419-e90ac5c316ab
# ╠═76c07283-c73b-4b3d-9396-e7f704bf7312
# ╠═99f4fcca-19bc-4a80-8472-198ed377abd0
# ╟─81c03358-1a48-4c33-802a-f2bfbf89f6d7
# ╟─0cbef7ea-c480-4679-88a2-f9e929f3c676
# ╠═b1b865d3-e4cd-4850-adba-9462b5930929
# ╟─0d55edfd-f003-466a-a12c-b4463f28168b
# ╠═bfb209de-c20b-455a-b45e-9f544d4877b0
# ╟─159b8d41-a5ce-494a-9566-df34547b1a8d
# ╟─7700b89e-8322-4d42-9b35-dd391731b5e3
# ╠═17dd25fd-1640-4f4b-9349-e35d73a897b8
# ╠═94e8d30c-9280-427c-a5bf-c3370712e69f
# ╠═d8f88426-d2a1-4ae1-8571-2ea80975a4fd
# ╟─e9c2610f-f94d-4426-be0f-d11d8e1498ab
# ╠═546844bf-3e36-4e9f-9bc9-aeb92649b108
# ╟─d628a011-a1fb-4947-8da1-2a1c7e2cf048
# ╠═4622c7a9-a35c-44ef-9c4a-e6d09faa6e08
# ╟─578f8606-4873-4c96-8112-1199fa486af1
# ╠═f0a5869c-e3eb-425d-aab4-02b8018ee6a3
# ╟─6ff92112-41d4-47bb-9eb0-4a2499f386fd
# ╠═c315d2b8-1a38-4ae1-b6d1-151b1a99cd20
# ╟─d4da9d8c-2064-4d1f-b397-4092ef0de660
# ╠═52a9641b-ea4c-45ce-b52f-ed3d47570a42
# ╟─3fc06973-839c-47d4-b0ca-2cdef798d5ca
# ╟─7c4e4208-944b-4c7e-a2fd-6275d2347a19
# ╟─fe1536be-8589-48ea-bbb8-9df739c4766c
# ╟─868a916b-1aa3-4f6d-acf8-84f9b955b3e5
# ╟─5c858f72-fcda-45b5-b9a3-036e5cc2d9e9
# ╠═1885527b-2fb8-46e1-99a5-a21e0efd0358
# ╟─cb3e7da3-4146-477e-8f24-6a55a1e800c7
# ╠═fdc4ba68-07b7-4dcf-8e6b-d4e55a4f8826
# ╠═664f2131-553d-4a89-ac5f-cf8138c59b4d
# ╟─c8bb9329-3c83-4cff-83bc-de79ab29c5a4
# ╠═f64da374-d99e-4be1-87ca-23261c56a26b
# ╠═30cbcfc4-0780-47d2-aec0-d990436e1693
# ╟─a1f98368-2ecb-40bc-80d2-b55d4224f2fc
# ╟─e0b1b009-e7a4-47dd-bd17-e32377b2d566
# ╟─13bdfa35-562b-49a5-978a-26c4a4920359
# ╠═f76d1f42-d5b9-40b5-9f31-667c04aee909
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

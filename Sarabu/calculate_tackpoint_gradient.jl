using LinearAlgebra, Flux

include(params.jl)

function tackpoint(x0, x_tar_c, maxiter, vel_wind, wind_dir, step_inc)
	#calculating initial xt---------------------------------------------------------
    vect_tar = (x_tar_c - x0) #generate vector to active target
    mid = 0.5 .* vect_tar + x0 #place point in the middle
    deviation = ((x_tar_c[2]-x0[2])^2/(x_tar_c[1]-x0[1]))/(atan(deg2rad(15))) #define deviation of point from vect_tar
    xt0 = [-deviation, 0] + mid #initial tack point
	
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
        pathtime_xtc = pathtime(x_tar_c, x0, xt, vel_wind, wind_dir, step_inc) # path time if current tack point is used
        pathtime_xtc_xplus = pathtime(x_tar_c, x0, xt + plength_x, vel_wind, wind_dir, step_inc)
        pathtime_xtc_xminus = pathtime(x_tar_c, x0, xt - plength_x, vel_wind, wind_dir, step_inc)
        pathtime_xtc_yplus = pathtime(x_tar_c, x0, xt + plength_y, vel_wind, wind_dir, step_inc)
        pathtime_xtc_yminus = pathtime(x_tar_c, x0, xt - plength_y, vel_wind, wind_dir, step_inc)

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
        minpathtime = pathtime(x_tar_c, x0, xt, vel_wind, wind_dir, step_inc) # calculate corresponding minpathtime
        minpathtime_p = pathtime(x_tar_c, x0, xt_p, vel_wind, wind_dir, step_inc)

        # prevent entering into infeasible region
        if minpathtime_p < minpathtime
            xt = xt_p
            minpathtime = minpathtime_p
            ax = 0.1 # reset acceleration
            ay = 0.1 # reset acceleration
        end
        println(iter)
        println(xt)
        println(minpathtime)
    end
    return xt
end

function tackpoint_GD(x0, x_tar_c, maxiter, vel_wind, wind_dir, step_inc)
	#calculating initial xt---------------------------------------------------------
    vect_tar = (x_tar_c - x0) #generate vector to active target
    mid = 0.5 .* vect_tar + x0 #place point in the middle
    deviation = ((x_tar_c[2]-x0[2])^2/(x_tar_c[1]-x0[1]))/(atan(deg2rad(15))) #define deviation of point from vect_tar
    xt0 = [-deviation, 0] + mid #initial tack point
	
	xt = xt0

    iter = 0

    while iter < maxiter # set maximum number of iterations
        # Define loss function
        loss(xt) = pathtime(x_tar_c, x0, xt, vel_wind, wind_dir, step_inc)

        # Compute gradients using Flux's gradient function
        grad_xt = Flux.gradient(x -> loss(x), xt)[1]

        # Update xt using gradient descent
        learning_rate = 0.1
        xt -= learning_rate * grad_xt

        iter += 1
		
        println(iter)
        println(xt)
    end
    return xt
end

xt_GD = tackpoint_GD(x0, x_tar_c, maxiter, vel_wind, wind_dir, step_inc)
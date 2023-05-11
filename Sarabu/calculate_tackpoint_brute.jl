using LinearAlgebra

my_acos(x) = x ≈ 1 ? zero(x) : x ≈ -1 ? one(x)*π : acos(x)

function pathtime(x_tar_c, x0, xt, vel_wind, wind_dir, step_inc)
    # x_tar_c = active target
    # x0 = current position
    # xt = tack coordinates

    theta_nogo = deg2rad(40) # deadzone
    velcons = 3 # velocity increase constant
    deg_int = 5 # constant provided by manufacture
    
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
        theta_rw = real(my_acos(dot(heading, wind_dir)/(norm(heading)*norm(wind_dir))))
        
        if dot(wind_dir, heading)/(norm(wind_dir)*norm(heading)) < 0
            # calculate velocity made good
            u = cos(theta_rd)*(vel_wind/cos(theta_nogo))*(1+0.01*velcons)^((abs(theta_rw-theta_nogo)*180)/(pi*deg_int))
            # check if within deadzone
            v = (rad2deg(my_acos(dot(heading,wind_dir)/(norm(heading)*norm(wind_dir))))<140) ? u : 0
        elseif dot(wind_dir,heading)/(norm(wind_dir)*norm(heading)) > 0
            # calculate velocity made good if downwind
            v = cos(theta_rd)*(vel_wind/cos(theta_nogo))*(1+0.01*velcons)^((abs(pi-theta_rw-theta_nogo)*180)/(pi*deg_int))
        end
        
        vel_tack = v/cos(theta_rd)
        
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
            v = (rad2deg(my_acos(dot(heading,wind_dir)/(norm(heading)*norm(wind_dir)))) < 140)*u
        elseif dot(wind_dir,heading)/(norm(wind_dir)*norm(heading)) > 0
            v = cos(theta_rd)*(vel_wind/cos(theta_nogo))*(1+0.01*velcons)^((abs(pi-theta_rw-theta_nogo)*180)/(pi*deg_int))
        end
        vel_tack = v/cos(theta_rd)
        if vel_tack > 10e-2
            time_increment_2 = 0.1/vel_tack
            time_elapsed_2 = time_elapsed_2 + time_increment_2
            xc = 0.1*heading/norm(heading) + xc
            step_c_2 = step_c_2 + 1
        elseif vel_tack < 10e-2
            time_elapsed_1 = 10e6
            break
        end
    end
    t_2 = time_elapsed_2
    t = t_1+t_2
    return t
end

function tackpoint(x0, x_tar_c, maxiter, vel_wind, wind_dir, start = [])
    #calculating initial xt---------------------------------------------------------
    vect_tar = (x_tar_c - x0) #generate vector to active target
    mid = 0.5 .* vect_tar + x0 #place point in the middle
    deviation = ((x_tar_c[2]-x0[2])^2/(x_tar_c[1]-x0[1]))/(atan(deg2rad(15))) #define deviation of point from vect_tar
    xt0 = [-deviation, 0] + mid #initial tack point

    if start == []
        xt = xt0
    else
        xt = start
    end

    iter = 0
    step_inc = 0.1 # step size

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

#Routine to determine the best tack point between two marks and calculate time
x_tar_c = [0.1,10] #active destination
x0 = [0,0] #starting point
maxiter = 100 #maximum iterations

vel_wind = 1 # wind velocity
wind_dir = [0, -1] # wind direction
start = [-5, 5]

xt = tackpoint(x0, x_tar_c, maxiter, vel_wind, wind_dir, start)
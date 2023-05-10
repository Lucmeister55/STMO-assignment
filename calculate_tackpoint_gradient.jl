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
    angle = my_acos(dot(normalize(vect_tar), normalize(wind_dir)))
    if angle < 45
        deviation = -0.1norm(vect_tar)
    else
        deviation = 0.1norm(vect_tar)
    end

    mid = 0.5 .* vect_tar + x0 #place point in the middle
    xt0 = [deviation, 0] + mid #initial tack point

    if start == []
        xt = xt0
    else
        xt = start
    end

    iter = 0
    step_size = 0.1 # step size

    while iter < maxiter # set maximum number of iterations
        # calculate path time at current tack point
        path_time = pathtime(x_tar_c, x0, xt, vel_wind, wind_dir, step_size)
        
        # approximate gradient using finite differences
        dx = [step_size, 0]
        dy = [0, step_size]
        grad_x = (pathtime(x_tar_c, x0, xt + dx, vel_wind, wind_dir, step_size) - path_time) / step_size
        grad_y = (pathtime(x_tar_c, x0, xt + dy, vel_wind, wind_dir, step_size) - path_time) / step_size
        
        # update tack point using gradient descent
        xt -= 0.01 * [grad_x, grad_y]
        
        iter += 1
        println(iter)
        println(xt)
        println(path_time)
    end
    
    return xt
end

#Routine to determine the best tack point between two marks and calculate time
x_tar_c = [0.1,10] #active destination
x0 = [0,0] #starting point
maxiter = 1000 #maximum iterations

vel_wind = 1 # wind velocity
wind_dir = [0, -1] # wind direction
start = [-5, 5]

xt = tackpoint(x0, x_tar_c, maxiter, vel_wind, wind_dir, start)
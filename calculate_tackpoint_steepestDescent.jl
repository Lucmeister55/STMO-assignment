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

function SteepestDescent_DLS(x)
    Xj=[0;0]
    Xs=x #Starting points
    Xi=[Xs[1];Xs[2]]
    Xt=Xi #Xt holds Xs information from previous iteration
    i=1 #Initialise iterations

    # function to compute the gradient at a given point (x,y)
    function grad(x,y)
        gx = -((1 + y^2)/x^2 - 2) * x^3
        gy = -2 * y / x^2 - 2 * x * y / (x^2 * (1 + y^2)^2)
        return [gx, gy]
    end

    # function to compute the objective function at a given point (x,y)
    function obj(x,y)
        return 12*x^2 + (x^2 + y^2 + 100)^2/(x^2 * y^2)
    end

    # convergence criteria to detect near zero gradient
    while abs.(grad(Xi[1],Xi[2])) .> [1e-7,1e-7]
        Xj=Xi
        gradi=grad(Xi[1],Xi[2]) #gradient of objective function at x,y
        al=10 #search bracket length
        a=0.5*al #initial search length to enter while statement
        ai=al
        il=0 #linesearch iteration
        as=0 #as is lower bracket of line search
        at=al #at is higher bracket of line search

        # convergence criteria to detect very small changes
        while abs(a-ai) > 1e-6
            ai=a
            # Dichotomous line search is employed with ap and aq used to
            # probe function value at either side of a length
            # the lengths are expressed as a function of the local Xj vector
            ap=a*0.99 #ap is probe length (at Xj)
            aq=a*1.01 #aq is probe length in opposite direction (Xj)
            # objp and objq are function values at either side of the search length
            objp=obj(Xj[1]-ap*gradi[1], Xj[2]-ap*gradi[2])
            objq=obj(Xj[1]-aq*gradi[1], Xj[2]-aq*gradi[2])
            # select side of the bracket with lower function value
            if objp<objq #if left side gives lower function value
                at=a #set right bracket as a
                as=as #leave left bracket as is
            elseif objp>objq #if right side gives lower function value
                at=at #leave right bracket as is
                as=a #set left bracket as a
            end
            a=0.5*(at-as) #set new search length as in between the brackets
            il=il+1 #update iteration for line search
        end

        Xj=Xi-a*gradi #update Xj vector
        Xt=Xi #update Xt
        Xi=Xj #update Xj
        i=i+1 #update iteration value
        #plot(Xj[1],Xj[2],'ro'); hold on #plot how solution travels
    end

    objs=obj(Xj[1],Xj[2]) #objective function at final solution
    Xk=[Xj;i;objs]
    ti=i #total iterations
    plot(i,Xj[1],'-'); 
    plot(i,Xj[2]);
    plot(i,obj(Xj[1],Xj[2]));
    m=minimum([0;minimum(x)]);
    n=maximum([10;maximum(x)]);
    m=collect(m:0.1:n);
    n=collect(m:0.1:n);
    M, N = meshgrid(m,n);
    O=(12+(M.^2)+(1+N.^2)./(M.^2)+(((M.*N).^2+100))./((M.*N).^4))*0.1;
    contour(M,N,O,[10e5,10e5,10e4,10e4,5000,5000,1000,1000,100,100,10,10,8,8,6,6,4,4,2,2]);
    title("Case 4 - Tack Point Optimization")
    xlabel("X")
    ylabel("Y")
end


#Routine to determine the best tack point between two marks and calculate time
x_tar_c = [0.1,10] #active destination
x0 = [0,0] #starting point
maxiter = 1000 #maximum iterations

vel_wind = 1 # wind velocity
wind_dir = [0, -1] # wind direction
start = [-5, 5]

xt = tackpoint(x0, x_tar_c, maxiter, vel_wind, wind_dir, start)
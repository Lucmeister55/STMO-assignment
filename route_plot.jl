using Dates
using Plots
using LinearAlgebra

function track(x0, x_tar, target_tol)
    xp = x0 # current position
    # Initialize target tracking
    num_targets = size(x_tar, 1) # get number of targets
    target_tracker = zeros(Int64, num_targets) # create file indicating no targets have been achieved
    itt = 1 # target tracker iteration
    while itt < num_targets + 1
        target_tracker[itt] = 0
        itt += 1
    end

    x_tar_c = [x_tar[1,1], x_tar[1,2]] # set initial target
    tar_act = 1 # active target is 1

    # Target tracking loop
    x = zeros(0)
    y = zeros(0)
    while norm(x_tar_c - xp) > target_tol / 10 # if target is not achieved
        step_direction = (x_tar_c - xp) / norm(x_tar_c - xp)
        stepsize = 0.1
        xp_p = xp # store previous location
        xp = step_direction .* stepsize .+ xp
        append!(x, xp[1])
        append!(y, xp[2])
        if norm(x_tar_c - xp) < target_tol
            if tar_act < size(x_tar, 1)
                target_tracker[tar_act] = 1
                x_tar_c = [x_tar[tar_act + 1, 1], x_tar[tar_act + 1, 2]]
                tar_act += 1
            end
        end
        println(norm(x_tar_c - xp))
    end
    return x, y
end

# plot single tack path
t0 = now() # record start time of function
x0 = [0,0] # starting points
target_tol = 4 # radius in metres within target to qualify
#x_tar = [-13.34947689 14.84067334; 20 50; -80.87069946 5.397111456; -100 0; -80.75569946 -7.316388104; 20 -50; 19.61212559 -49.87259138; 0 0] # specify targets in io
#xt = [-13.34947689 14.84067334]
xt = [-3.85 4.130434782608695]
x_tar = [xt; 0.1 10]
wind_dir = [0, -1]

x, y = track(x0, x_tar, target_tol)

scatter(x, y)
# Plot the wind vector
xw = [x_tar[2, 1], x_tar[2, 1]+wind_dir[1]]
yw = [x_tar[2, 2], x_tar[2, 2]+wind_dir[2]]
plot!(xw,yw, arrow=true, arrowsize=0.5)
savefig("trajectory.png")

# Record elapsed time
solvertime = (now() - t0)
println(solvertime)

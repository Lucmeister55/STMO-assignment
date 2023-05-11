using Distributions

# Define the parameters
d = 2   # Dimension of the continuous component of the state
NI = 2  # Number of possible discrete states
U = [-1, 0, 1]  # Set of controls
Tf = 10  # Final time
dt = 0.01  # Time step size

# Define the continuous dynamics function
function f(x, q, u)
    if q == 1
        return [x[1] + u[1]*dt, x[2] + u[2]*dt]
    else
        return [x[1] + u[2]*dt, x[2] + u[1]*dt]
    end
end

# Define the switching times
t1 = 3
t2 = 6
t3 = 9

# Define the piecewise constant function Q
function Q(t)
    if t < t1
        return 1
    elseif t < t2
        return 2
    elseif t < t3
        return 1
    else
        return 2
    end
end

# Define the initial conditions
x0 = [0.0, 0.0]
q0 = 1

# Define the standard Brownian motion process
W = zeros(d, Int(Tf/dt))
dW = sqrt(dt)*randn(d, Int(Tf/dt))
for i = 2:Int(Tf/dt)
    W[:,i] = W[:,i-1] + dW[:,i-1]
end

# Simulate the system
x = zeros(d, Int(Tf/dt))
q = zeros(Int(Tf/dt))
x[:,1] = x0
q[1] = q0
for i = 2:Int(Tf/dt)
    u = rand(U)
    x[:,i] = f(x[:,i-1], q[i-1], u) + dW[:,i-1]
    q[i] = Q(i*dt)
end

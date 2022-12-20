#Implement the 3 body problem dynamics
using ForwardDiff
using LinearAlgebra

μ = 1.215e-2 #gravitational parameter
μ1 = 1-μ
μ2 = μ

pose_m1 = [-μ, 0, 0]
pose_m2 = [1-μ, 0, 0]

L = 3.850e5 #in km - distance between centers of m1 and m2
V = 1.025 #in km/s - orbital velocity of m1
T = 2.361e6 #in seconds - orbital period of m1 and m2

h_s = 100 #time step for RK4 in seconds # original 10 seconds

time_scale = T/(2*pi)

h_rk4 = h_s/time_scale #non dimensionalized time

#Determine amount of knot points
one_min = 60/h_s #number of timesteps for 1 minute
min_day = 1440 #number of minutes this is a day
num_days = 3 #120 # 3 days original
t_f=num_days*min_day 
Nt = Int(one_min*t_f) #time steps for duration


function potential_energy(X)
    
    x = X[1]
    y = X[2]
    z = X[3]
   
   
    r1 = sqrt((x+μ2)^2 + y^2 + z^2) 
    r2 = sqrt((x-μ1)^2 + y^2 + z^2)
    #assuming m3 is unit mass
   
    U = (-μ1/r1)-(μ2/r2)-0.5*(μ1*μ2)
   
    return U
end

function effective_potential(X)
    
    x = X[1]
    y = X[2]
    z = X[3]
   
   
    #r1 = sqrt((x+μ2)^2 + y^2) 
    #r2 = sqrt((x-μ1)^2 + y^2)

    r1 = sqrt((x+μ2)^2 + y^2 + z^2) 
    r2 = sqrt((x-μ1)^2 + y^2 + z^2)
    #assuming m3 is unit mass
   
   U = (-μ1/r1)-(μ2/r2)-0.5*(x^2+y^2)
   
    return U
end

function CR3BPdynamics(rv) #Three body dynamics in Sun-Earth System

    #changed for EARTH MOON !!!
    μ = 1.215e-2    
    r₁³= ((rv[1] + μ)^2.0     + rv[2]^2.0 + rv[3]^2.0)^1.5; # distance to m1, LARGER MASS
    r₂³= ((rv[1] - 1 + μ)^2.0 + rv[2]^2.0 + rv[3]^2.0)^1.5; # distance to m2, smaller mass
    # r₁³= ((x + μ)^2     + y^2 + z^2)^1.5; # distance to m1, LARGER MASS
    # r₂³= ((x - 1 + μ)^2 + y^2 + z^2)^1.5; # distance to m2, smaller mass

#     rvdot = zeros(6)
    rvdot = zeros(eltype(rv),6)
    rvdot[1:3] = [rv[4];rv[5];rv[6]]
    rvdot[4] = -((1.0 - μ)*(rv[1] + μ)/r₁³) - (μ*(rv[1] - 1.0 + μ)/r₂³) + 2.0*rv[5] + rv[1];
    rvdot[5] = -((1.0 - μ)*rv[2]      /r₁³) - (μ*rv[2]          /r₂³) - 2.0*rv[4] + rv[2];
    rvdot[6] = -((1.0 - μ)*rv[3]      /r₁³) - (μ*rv[3]          /r₂³);
    return rvdot
end



# function threebp_dynamics(x)
        
#     q = x[1:3]
#     v = x[4:6]
#     #a = zeros(3)
#     ẋ = zeros(eltype(x),6)
#     U_q = zeros(eltype(x),3)
#     U_q = (ForwardDiff.gradient(_x -> effective_potential(_x), q))
        
#     a[1] = 2*v[2] - U_q[1]
#     a[2] = -2*v[1] - U_q[2]
#     a[3] = -U_q[3]
    
#     ẋ = [v; a] #x dot is velocity and acceleration
    
#     return ẋ
# end

function RK4_satellite_potential(x, h)

    f1 = threebp_dynamics(x)
    f2 = threebp_dynamics(x+0.5*h_rk4*f1)
    f3 = threebp_dynamics(x+0.5*h_rk4*f2)
    f4 = threebp_dynamics(x+h_rk4*f3)
    
    #forward in time
    xnext = x+(h_rk4/6.0)*(f1+2*f2+2*f3+f4)
    
    #backward in time
    #xnext = x-(h_rk4/6.0)*(f1+2*f2+2*f3+f4)
    
    return xnext
    
end

function threebp_dynamics(x)
        
    q = zeros(eltype(x),3)
    v = zeros(eltype(x),3)
    q = x[1:3]
    v = x[4:6]
    #a = zeros(3)

    ẋ = zeros(eltype(x),6)
    U_q = zeros(eltype(x),3)
    
    U_q = (ForwardDiff.gradient(_x -> effective_potential(_x), q))
        
    #a[1] = 2*v[2] - U_q[1]
    #a[2] = -2*v[1] - U_q[2]
    #a[3] = -U_q[3]

    ẋ[1:3] = v
    ẋ[4] = 2*v[2] - U_q[1]
    ẋ[5] = -2*v[1] - U_q[2]
    ẋ[6] = -U_q[3]
    
    #ẋ = [v; a] #x dot is velocity and acceleration
    
    return ẋ
end



function get_reference_traj(x0, h)

    xk = zeros(6, Nt)

    xk[:,1] = x0

    for k=1:Nt-1
    
        xk[:,k+1] = RK4_satellite_potential(xk[:,k], h)
        
    end

    return xk

end


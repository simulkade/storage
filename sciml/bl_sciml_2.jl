# sciml practice on BL
# my first attempt at using the sciml infrastructure
# of Julia to do parameter estimation
using DifferentialEquations, Plots, Flux, Optim, DiffEqFlux #, DiffEqSensitivity

# define the domain (global variables)
L = 0.1 # meter
n = 100 # number of cells
dx = L/n
μ_w = 1e-3 # Pa.s
μ_o = 1.2e-3 # Pa.s
krw0 = 0.3
kro0 = 0.8
nw = 4.0
no = 2.0
swc = 0.1
sor = 0.2
sw_in = 1.0
sw_init = swc
u = 1e-5 # m/s
φ = 0.4 # porosity
k = 0.01e-12 # m^2 permeability

"""
water relative permeability
"""
function krw(sw, krw0, nw, swc, sor)
    if sw<swc
        res=0.0
    elseif sw>1-sor
        res = krw0
    else
        res = krw0*((sw-swc)/(1-sor-swc))^nw
    end
    return res
end # krw

"""
Oil relative permeability
"""
function kro(sw, kro0, no, swc, sor)
    if sw<swc
        res = kro0
    elseif sw>1-sor
        res = 0.0
    else
        res = kro0*((1-sw-sor)/(1-swc-sor))^no
    end
    return res
end #kro


"""
BL_ode converts the Buckley-Leverett equation to 
a set of ODEs using method of lines
dSW/dt + u/φ*dfw/dx = 0
p = [krw0, kro0, nw, no, swc, sor, muw, muo]
"""
function BL_ode(dsw, sw, p, t)
    krw0, kro0, nw, no, swc, sor = p
    # upwind
    sw_face = [sw_in; sw]

    # central
    # sw_face = [sw_in; 0.5*(sw[1:end-1]+sw[2:end]); sw[end]]

    λ_w = krw.(sw_face, krw0, nw, swc, sor)./μ_w
    λ_o = kro.(sw_face, kro0, no, swc, sor)./μ_o 
    fw = λ_w./(λ_w.+λ_o)
    
    for i in 1:length(sw)
        dsw[i] = -u/φ*(fw[i+1]-fw[i])/dx
    end
end # BL_ode()

function recovery_factor(sol)
    oil_init = 1.0 - dx/L*sum(sw0)
    R = (oil_init .- (1.0 .- dx/L.*sum(Array(sol), dims=1))')./oil_init
    return R[:]
end # recovery_factor

function pressure_drop(sol)
    λ_w = krw.(Array(sol), krw0, nw, swc, sor)./μ_w
    λ_o = kro.(Array(sol), kro0, no, swc, sor)./μ_o
    dp = dx*sum(u./(k.*(λ_w .+ λ_o)), dims=1)
    return dp[:]
end


# initial conditions
sw0 = sw_init.*ones(n)

# parameters
# krw0, kro0, nw, no, swc, sor, μ_w, μ_o, sw_in, u, φ, dx 
p0 = [krw0, kro0, nw, no, swc, sor]
tspan = (0.0, 6*3600)

prob = ODEProblem(BL_ode, sw0, tspan, p0)
sol = solve(prob)

plot(sol)

# extract the recovery factors
oil_init = 1.0 - dx/L*sum(sw0)
R = zeros(size(sol.t))
for i in 1:length(sol.t)
    R[i] = (oil_init-(1.0 - dx/L*sum(sol.u[i])))/oil_init
end

# calculate pressure drop:
# -(λ_w + λ_o)∇p = u
# -(λ_w + λ_o)dp/dx = u
dp = zeros(size(sol.t))
for i in 1:length(sol.t)
    λ_w = krw.(sol.u[i], krw0, nw, swc, sor)./μ_w
    λ_o = kro.(sol.u[i], kro0, no, swc, sor)./μ_o
    dp[i] = dx*sum(u./(k.*(λ_w .+ λ_o)))
end

# synthetic experimental data:
R_exp = R[1:floor(Int, length(sol.t)/10):length(sol.t)]
dp_exp = dp[1:floor(Int, length(sol.t)/10):length(sol.t)]
dp_max_exp = maximum(dp_exp)
t_exp = sol.t[1:floor(Int, length(sol.t)/10):length(sol.t)]
# scatter(t_exp, R_exp)

function loss(p)
    sol = solve(prob, Tsit5(), p=p, saveat = t_exp)
    # Note: the following loop is not differentiable and results in an error from Zygote
    # R = zeros(size(sol.t))
    # for i in 1:length(sol.t)
    #     loss += (R_exp[i] - (oil_init-(1.0 - dx/L*sum(sol.u[i])))/oil_init)^2
    # end

    # Note: the following form is differentiable:
    # loss = sum(abs2, (oil_init .- (1.0 .- dx/L.*sum(sol, dims=1))')./oil_init .- R_exp)
    loss_recovery = sum(abs2, recovery_factor(sol).-R_exp)
    loss_dp       = sum(abs2, (pressure_drop(sol).-dp_exp)./dp_max_exp)
    loss_total = loss_recovery+loss_dp
    
    return loss_total, sol
end
  
callback = function (p, l, pred)
    display(l)
    # plt = plot(pred, ylim = (0, 6))
    # display(plt)
    # Tell sciml_train to not halt the optimization. If return true, then
    # optimization stops.
    return false
end

# [krw0, kro0, nw, no, swc, sor]
p_est = [0.8, 0.3, 2.0, 2.0, 0.1, 0.15]

lb = [0.05, 0.05, 1.0, 1.0, 0.02, 0.02]
ub = [1.0, 1.0, 4.0, 4.0, 0.3, 0.3]      

function fit_model()
result_ode = DiffEqFlux.sciml_train(loss, p_est,
                                    ADAM(0.1),
                                    cb = callback,
                                    maxiters = 100,
                                    lower_bounds=lb,
                                    upper_bounds=ub)
    return result_ode
end                      

# p_min = [0.6875461937595114,
#  0.48905916428640167,
#  2.4236909544499854,
#  2.2117465258145685,
#  0.4544348486030391,
#  0.19496540359779663]

result_ode = fit_model()
p_min = result_ode.minimizer

prob = ODEProblem(BL_ode, sw0, tspan, p_min)
sol = solve(prob)
R = zeros(size(sol.t))
for i in 1:length(sol.t)
    R[i] = (oil_init-(1.0 - dx/L*sum(sol.u[i])))/oil_init
end

scatter(t_exp, R_exp)

plot!(sol.t, R)

scatter(t_exp, dp_exp)
plot!(sol.t, pressure_drop(sol))

using JFVM, Zygote, JFVMvis, ForwardDiff

# diffusion and reaction:
# dc/dt - D d2c/dx2 -kc = 0
# left BC: c=1
# right BC: dc/dx = 0
n_steps = 10
dt = 1000 # s
t_final = n_steps*dt

D = 1e-5 # m^2/s
k = 1e-4 # 
c0 = 0.0

m = createMesh1D(10, 1.0)

BC = createBC(m)

BC.left.a .= 0.0
BC.left.b .= 1.0
BC.left.c .= 1.0
Mbc, RHSbc = boundaryConditionTerm(BC)

D_cell = createCellVariable(m, D)
D_face = harmonicMean(D_cell)
Md = diffusionTerm(D_face)


c_init = createCellVariable(m, c0, BC)

function diff_react(k_reac)
    k_cell = createCellVariable(m, k_reac[1])
    Ms = linearSourceTerm(k_cell)
    for t in dt:dt:t_final
        Mt, RHSt = transientTerm(c_init, dt)
        # c = solveLinearPDE(m, Mt+Mbc-Md-Ms, RHSt+RHSbc)
        c = (Mt+Mbc-Md-Ms)\(RHSt+RHSbc)
        c_init.value .= c
    end
    return sum(c_init.value) # return a scalar to test AD methods
end

# gradient(diff_react, k) # does not work :-(

ForwardDiff.gradient(diff_react, [k])
# visualizeCells(c_init)

# conclusion:
# Zygote does not work with my JFVM package
# perhaps it needs a lot of changes that is not currenlty affordable

# ForwardDiff goes to some stage but gets stock with a linear solver error

# It seems that the only way is to implement things from scratch, which is
# going to take time
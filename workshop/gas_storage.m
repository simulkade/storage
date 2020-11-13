% simple gas storage example
% density data for H2 and methane
clc
T_ref = 70+273.15; % K reservoir temperature
p_ref = 200e5; % Pa reservoir pressure
rho_ref = py.CoolProp.CoolProp.PropsSI('D','P',p_ref,'T',T_ref,'CH4');
mu_ref = py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',p_ref,'T',T_ref,'CH4');

n_data = 20;
p_range = linspace(p_ref/2, p_ref*1.5, n_data)';
rho_data = zeros(n_data, 1);
comp = 'CH4';
for i=1:n_data
    rho_data(i) = py.CoolProp.CoolProp.PropsSI('D','P',p_range(i),'T',T_ref, comp);
end

rho_fit = polyfit(p_range, rho_data, 1); % line for now
drho_fit = polyder(rho_fit);
rho = @(p)polyval(rho_fit, p);
drho = @(p)polyval(drho_fit, p);

figure(1);
plot(p_range, rho_data, 'o', p_range, polyval(rho_fit,p_range));
close all;


% set up the simulation
perm = 0.01e-12; % m2
poros = 0.4; % chalk porosity
R_field = 500; % m field radius
H_field = 30; % m field thickness
well_radius = 0.05; % m
Nx = 100;
Ny = 20;
m = createMeshCylindrical2D(Nx, Ny, R_field, H_field);
% shift the mesh to create the well
m.cellcenters.x = m.cellcenters.x+well_radius;
m.facecenters.x = m.facecenters.x+well_radius;


% gas production
surplus_elec = 100e6; % Watt
eta_transmission = 0.95; % conversion and transmission efficiency for electricity
eta_electrolysis = 0.7; % electrolysis efficiency
H2_exergy = 200000; % J/mol H2
H2_prod = surplus_elec*eta_transmission*eta_electrolysis/H2_exergy; % mol/s H2
CH4_prod = H2_prod/2; % assuming 100% efficiecny for H2 to CH4 and free carbon source
MW_CH4 = 0.016; % kg/mol
CH4_flow = CH4_prod*MW_CH4/rho_ref; % m3/s CH4 flow at the reference reservoir pressure and temperature (reservoir)
well_area = pi*well_radius*2*H_field; % m2 well area
u_inj = CH4_flow/well_area; % m/s injection velocity

BC = createBC(m);
BC.left.a(:) = perm/mu_ref; BC.left.b(:)=0.0; BC.left.c(:) = -u_inj;
BC.right.a(:) = 0.0; BC.right.b(:)=1.0; BC.right.c(:) = p_ref;
[M_BCp, RHS_BCp] = boundaryCondition(BC);

% Tracer
D_val = 1e-6; % m2/s diffusion coefficient
D_cell = createCellVariable(m, D_val);
D_face = harmonicMean(D_cell);
M_dc = diffusionTerm(D_face);

BCc = createBC(m);
BCc.left.a(:) = 0.0; BCc.left.b(:)=1.0; BCc.left.c(:)=1.0;
[M_bcc, RHS_bcc] = boundaryCondition(BCc);

c0 = 0.0;
c_init = createCellVariable(m, c0);

V_dp = 0.2; % Dykstra-Parsons coeff
clx = 1.1; % correlation length x direction
cly = 0.1; % correlation length y direction
rng(1);
perm_field = field2d(Nx, Ny, perm, V_dp, clx, cly);
k = createCellVariable(m, perm_field);
k_face = harmonicMean(k);
mobil = k_face./mu_ref;
hFig = figure(1);
set(hFig, 'Position', [200 200 1500 500]);
visualizeCells(k); axis normal; shading interp;

% initial condition
p_init = createCellVariable(m, p_ref, BC); % ignoring gravity
p_val = p_init;
rho_val = createCellVariable(m, rho_ref);
drho_val = createCellVariable(m, 0);
rho_init = rho_val;

% injecting for a week
t_final = 7*24*3600; % s injection time
dt = 3600; % s time step

t=0;
while t<t_final
    t=t+dt;
    for i = 1:3 % internal loop
        % solve pressure equation
        drho_val.value = drho(p_val.value);
        [M_trans, RHS_trans] = transientTerm(p_init, dt, poros.*drho_val);
        %  RHS_ddt_p   = constantSourceTerm(poros/dt*(rho_val - rho_init));
        p_face = arithmeticMean(p_val);
        grad_p = gradientTerm(p_val);
        drho_face = funceval(drho, p_face);
        rho_face = funceval(rho, p_face);

        RHS_source = divergenceTerm(-mobil.*drho_face.*p_face.*grad_p);
        M_conv = convectionTerm(-mobil.*drho_face.*grad_p);
        M_diff = diffusionTerm(-mobil.*rho_face);
        % rho_face = funceval(@(x)polyval(rho_fit, x), p_face); % arithmeticMean(rho_val);
        % M_diff_p    = diffusionTerm(rho_face.*mobil);
        p_val       = solvePDE(m, M_trans+M_diff+M_conv+M_BCp, RHS_BCp+RHS_trans+RHS_source);

        % velocity vector
        

    %         # solve heat equation
    %         α                  = poros*rho_val*cp_water+(1-poros)*
    %             rho_rock*cp_rock
    %         M_trans, RHS_trans = transientTerm(theta_init, dt, α)
    %         M_conv             = convectionTerm(cp_water*rho_face.*u)
    %         theta_val          = solveLinearPDE(m, 
    %             M_BCt + M_conv + M_trans + M_conductivity,
    %             RHS_BCt + RHS_trans)
    %         
        % update density values
        % rho_val.value(:) = polyval(rho_fit, p_val.value(:));
    end % end of inner loop
    p_init = p_val;
    %visualizeCells(p_val);drawnow;
    
    u = -mobil.*gradientTerm(p_val);
    % solve for concentration
    [M_tc, RHS_tc] = transientTerm(c_init, dt, poros);
    M_cc = convectionUpwindTerm(u);
    c_init = solvePDE(m, M_tc+M_cc+M_bcc-M_dc, RHS_tc+RHS_bcc);
    
    visualizeCells(c_init);axis normal;drawnow;
end

% visualizeCells(p_val);
hFig = figure(2);
set(hFig, 'Position', [200 200 1500 500]);
visualizeCells(p_val); shading interp; axis normal;

hFig = figure(3);
set(hFig, 'Position', [200 200 1500 500]);
visualizeCells(c_init); shading interp; axis normal;

% copy-paste from my julia code
% for i in 1:3 # internal loop
%         # solve pressure equation
%         RHS_ddt_p   = constantSourceTerm(poros/dt*(rho_val - rho_init))
%         water_mobil = harmonicMean(perm_field./mu_val)
%         rho_face = arithmeticMean(rho_val)
%         M_diff_p    = diffusionTerm(-rho_face.*water_mobil)
%         p_val       = solveLinearPDE(m, M_diff_p + M_BCp,
%             RHS_BCp - RHS_ddt_p)
%         
%         # velocity vector
%         u = -water_mobil.*gradientTerm(p_val)
%         
%         # solve heat equation
%         α                  = poros*rho_val*cp_water+(1-poros)*
%             rho_rock*cp_rock
%         M_trans, RHS_trans = transientTerm(theta_init, dt, α)
%         M_conv             = convectionTerm(cp_water*rho_face.*u)
%         theta_val          = solveLinearPDE(m, 
%             M_BCt + M_conv + M_trans + M_conductivity,
%             RHS_BCt + RHS_trans)
%         
%         # update density and viscosity values
%         rho_val.value[:] = polyval(rho_fit, theta_val.value + T0)
%         mu_val.value[:]  = polyval(mu_fit, theta_val.value + T0)
%     end # end of inner loop

% simple gas storage example
% density data for H2 and methane
clc

% extract injection data:
date_begin = '2020-10-01';
date_end = '2021-01-20';
d_long = extract_elect_data(date_begin, date_end);
d_long.HourDK = datetime(d_long.HourDK, 'Format', 'yyyy-MM-dd''T''HH:mm:ss');
date_num = datenum(d_long.HourDK);
future_demand = 1.5; % including 50% increase of demand for heating and cars
future_capacity = 4.0; % windmill capacity will be quadruppled
[wind_elec, elec_demand, surplus, shortage] = extract_surplus(d_long, future_capacity, future_demand);

% Lets inject into the reservoir for 4 month and produce/inject later for 2 months
end_storage = '2020-11-30';
ind_end_storage = find(datenum(end_storage)<date_num, 1); % index of end of storage

% injecting for a week
t0 = date_num(1)*24*3600; % [s]
t = t0;
t_inj = date_num(ind_end_storage)*24*3600; % s injection time
% t_prod = 7*24*3600; % s production time
t_final = date_num(end)*24*3600; % total simulation time
dt = 3600; % s time step (one hour time steps)


% Input data
% set up the simulation
perm = 0.01e-12; % m2
poros = 0.18; % sandstone porosity
R_field = 1000; % m field radius
H_field = 100; % m field thickness
well_radius = 5*0.0254/2; % m (converting inch to m)
frac_aperture = 0.01; % 1 mm fracture (I know too big :-))
N_well = 1; % number of wells
Nx = 60;
Ny = 20;
x_mesh = linspace(0, R_field, Nx);
y_mesh1 = linspace(0, H_field, Ny);
y_mesh2 = y_mesh1+frac_aperture;
y_mesh = zeros(1,2*Ny);
y_mesh(1:2:2*Ny-1) = y_mesh1;
y_mesh(2:2:2*Ny) = y_mesh2;

p_0 = 75e5; % [Pa] for a depleted reservoir
p_std = 1e5; % [Pa]

eta_comp = 0.8; % compressor efficiency
eta_turbine = 0.8; % turbine efficiency
eta_gen = 0.9; % generator efficiency

% stimulation radius
r_stim = 100*well_radius; % stimulation radius
k_stim = perm; % stimulated permeability

comp = 'Air';
T_ref = 135+273.15; % K reservoir temperature
p_ref = 3500e5; % Pa reservoir pressure
rho_ref = py.CoolProp.CoolProp.PropsSI('D','P',p_ref,'T',T_ref,comp); % kg/m3
mu_ref = py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',p_ref,'T',T_ref,comp); % Pa.s

n_data = 100;
p_range = linspace(p_std, p_ref*3, n_data)';
rho_data = zeros(n_data, 1);
for i=1:n_data
    rho_data(i) = py.CoolProp.CoolProp.PropsSI('D','P',p_range(i),'T',T_ref, comp);
end

rho = @(p)interp1(p_range, rho_data, p); % line for now
% drho_fit = polyder(rho_fit);
% rho = @(p)polyval(rho_fit, p);
eps1 = 10; % [Pa]
drho = @(p)((rho(p+eps1)-rho(p))/eps1);

figure(1);
plot(p_range, rho_data, 'o', p_range, rho(p_range));
close all;



m = createMeshCylindrical2D(x_mesh, y_mesh);
% shift the mesh to create the well
m.cellcenters.x = m.cellcenters.x+well_radius;
m.facecenters.x = m.facecenters.x+well_radius;


% % gas production
% surplus_elec = 10e6; % 10e6 Watt (10 MW)
eta_transmission = 0.95; % conversion and transmission efficiency for electricity

well_area = pi*well_radius*2*H_field; % m2 well area


BC = createBC(m);

BC.right.a(:) = 0.0; BC.right.b(:)=1.0; BC.right.c(:) = p_0;
% BC.right.a(:) = 1.0; BC.right.b(:)=0.0; BC.right.c(:) = 0;

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
perm_field = field2d(Nx-1, 2*Ny-1, perm, V_dp, clx, cly);
% perm_field(1:2,:) = k_stim; % stimulation
perm_field(1:Nx/2, 3:4:2*Ny-1) = 100*perm; % fracture permeability
k = createCellVariable(m, perm_field);
k_face = harmonicMean(k);
mobil = k_face./mu_ref;
hFig = figure(1);
set(hFig, 'Position', [200 200 1500 500]);
visualizeCells(k); % axis normal; shading interp;

% initial condition
p_init = createCellVariable(m, p_0, BC); % ignoring gravity
p_val = p_init;
rho_val = createCellVariable(m, rho_ref);
drho_val = createCellVariable(m, 0);
rho_init = rho_val;


p_hist = zeros(length(date_num), 1); % array to store average injection pressure
t_hist = zeros(length(date_num), 1);
m_hist = zeros(length(date_num), 1);
j=0;
T_inj = 135+273.15; % K injection temperature
p_atm = 1.0e5;
p_inj = p_ref; % initial estimate
MW_air = 0.029; % kg/mol
comp_ratio = 3.5;
elec_max = 10e6; % J/s
shortage_max = elec_max; % J/s
while t<t_final
    t=t+dt;
    j = j+1;
    if t<t_inj % time dependent injection boundary condition
        w_isentropic = compressor_station(p_atm, p_inj, comp_ratio, comp, T_inj);
        w_real = w_isentropic/eta_comp;
        rho_inj = py.CoolProp.CoolProp.PropsSI('D','P',p_inj,'T',T_inj, comp); % density kg/m3, must change with pressure
        m_inj = (surplus(j)>0)*elec_max*eta_transmission/w_real*MW_air; % kg/s
        u_inj = m_inj/rho_inj/well_area/N_well; % injection velocity
        BC.left.a(:) = perm/mu_ref; BC.left.b(:)=0.0; BC.left.c(:) = -u_inj;
        [M_BCp, RHS_BCp] = boundaryCondition(BC);
    elseif (t>t_inj)
        % production/storage time
%         BC.left.a(:) = perm/mu_ref; BC.left.b(:)=0.0; BC.left.c(:) = u_inj;
%         [M_BCp, RHS_BCp] = boundaryCondition(BC);
        
        w_isentropic = compressor_station(p_atm, p_inj, comp_ratio, comp, T_inj);
        w_real = w_isentropic/eta_comp;
        rho_inj = py.CoolProp.CoolProp.PropsSI('D','P',p_inj,'T',T_inj, comp); % density kg/m3, must change with pressure
        if surplus(j)>0
            elec_bound = elec_max;
        else
            elec_bound = shortage_max;
        end
        m_inj = sign(surplus(j)-shortage(j))*elec_bound*eta_transmission/w_real*MW_air; % kg/s
        u_inj = m_inj/rho_inj/well_area/N_well; % injection velocity
        BC.left.a(:) = perm/mu_ref; BC.left.b(:)=0.0; BC.left.c(:) = -u_inj;
        [M_BCp, RHS_BCp] = boundaryCondition(BC);
    end
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
    
    
    p_hist(j) = mean(p_face.xvalue(1,:));
    t_hist(j) = t;
    p_inj = p_hist(j);
    m_hist(j) = m_inj;
    disp(p_inj);
    
    u = -mobil.*gradientTerm(p_val);
    % solve for concentration
    [M_tc, RHS_tc] = transientTerm(c_init, dt, poros);
    M_cc = convectionUpwindTerm(u);
    c_init = solvePDE(m, M_tc+M_cc+M_bcc-M_dc, RHS_tc+RHS_bcc);
    
    % visualizeCells(c_init);axis normal;drawnow;
end

% visualizeCells(p_val);
hFig = figure(2);
set(hFig, 'Position', [200 200 1500 500]);
visualizeCells(p_val); shading interp; axis normal;

hFig = figure(3);
set(hFig, 'Position', [200 200 1500 500]);
visualizeCells(c_init); shading interp; axis normal;

figure(4);
plot(t_hist/(3600*24), p_hist/1e5);
xlabel('time[day]')
ylabel('bar')

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

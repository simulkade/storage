% testing the correlation \rho = \rho_0 (c_f \exp(p-p_0)) for different fluids
clc

%% Step 1: test the correlation for a liquid
component = 'NH3';
T_ref = 70+273.15; % [K] temperature in the North Sea reservoirs
p0 = 100e5; % [Pa] pressure of a depleted reservoir
p1 = 300e5; % [Pa] pressure of a pressurized reservoir
p_range = linspace(p0, p1, 50); % pressure range

rho_eos = arrayfun(@(x)density(x, T_ref, component), p_range); % density data from coolprop
rho_cor = density_correlation(p0, p1, T_ref, component); % density-pressure correlation

rho_est = arrayfun(rho_cor, p_range);

plot(p_range, rho_eos, 'o', p_range, rho_est)
xlabel('Pressure [Pa]')
ylabel([component ' density [kg/m^3]'])
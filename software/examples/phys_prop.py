from CoolProp.CoolProp import PropsSI
import numpy as np
import pandas as pd


comp = 'Air'
T_ref = 135+273.15 # K reservoir temperature
p_ref = 3500e5     # Pa reservoir pressure
p_std = 1e5        # [Pa]
rho_ref = PropsSI('D','P',p_ref,'T',T_ref,comp)          # kg/m3
mu_ref  = PropsSI('VISCOSITY','P',p_ref,'T',T_ref,comp)  # Pa.s

n_data = 100
p_range = np.linspace(p_std, p_ref*3, n_data)
# rho_data = np.zeros(n_data)
rho_data = PropsSI('D','P',p_range,'T',T_ref, comp)
data_out = np.array([p_range, rho_data])
data_out = data_out.T
df = pd.DataFrame(data_out, columns=['pressure', 'density'])
df.to_csv('p_dens.csv')
# data_out.tofile('p_dens.csv', sep=',')
# df.to_csv
# for i in 1:n_data
#     rho_data(i) = py.CoolProp.CoolProp.PropsSI('D','P',p_range(i),'T',T_ref, comp);
# end

# rho = @(p)interp1(p_range, rho_data, p); % line for now
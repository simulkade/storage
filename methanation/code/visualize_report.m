% visualizing the methanation results
load methanation_zvi.mat
zvi_res = reactants;
dt = 10000; %[s]
t_end = 1000*dt; %[s]
t = (dt:dt:t_end)/(3600*24);

load methanation_natural.mat
meth_res = reactants;

hold all
h = figure(1)
yyaxis left
plot(t, zvi_res(:,1), '-', 'linewidth', 2)
plot(t, meth_res(:,1), '--', 'linewidth', 2)
ylabel('oil [mol]')


yyaxis right
plot(t, zvi_res(:,4), '-r', 'linewidth', 3)
plot(t, meth_res(:,4), '--r', 'linewidth', 3)
plot(t, zvi_res(:, 7), 'r-');
plot(t, meth_res(:, 7), '--r');
plot(t, 500*zvi_res(:,8), '-', 'color', 'magenta') 
plot(t, 500*meth_res(:,8), '--', 'color', 'magenta')
ylabel('CO2, CH4 [mol]')

ax = gca
ax.XAxis.FontSize = 14;
ax.YAxis(1).FontSize = 14;
ax.YAxis(2).FontSize = 14;
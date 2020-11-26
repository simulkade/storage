% run example 11 of phreeqc in PhreeqcMatlab
[phreeqc_rm, c_x_t] = PhreeqcAdvection('monod.phr' , 'methanation_1d.pqm');
plot(squeeze(c_x_t(:, 4, end))'); % plot the concentrations profile
comps = phreeqc_rm.GetComponents();
legend(comps{4:end});
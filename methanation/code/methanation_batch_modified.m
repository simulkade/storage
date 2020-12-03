% run example 11 of phreeqc in PhreeqcMatlab
phreeqc_rm = PhreeqcSingleCell('biomass_appelo_batch_modified.phr' , 'llnl.dat');

h_out = phreeqc_rm.GetSelectedOutputHeadings(1);

dt = 3600; %[s]
t_end = 100*dt; %[s]
t = dt:dt:t_end;
n_data = length(dt:dt:t_end);
reactants = zeros(n_data, length(h_out));

for i = 1:length(t)
    phreeqc_rm.RM_SetTime(t(i));
    phreeqc_rm.RM_SetTimeStep(dt);

    status = phreeqc_rm.RM_RunCells();
    % c_out = phreeqc_rm.GetConcentrations();

    reactants(i, :) = phreeqc_rm.GetSelectedOutput(1);
    
end
figure(1)
yyaxis left
plot(t, reactants(:,1))
yyaxis right
plot(t, reactants(:,2))

figure(2)
yyaxis left
plot(t, reactants(:, 4));
yyaxis right
plot(t, reactants(:, 7));

figure(3)
plot(t, reactants(:, 6));
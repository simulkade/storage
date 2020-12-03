% run example 11 of phreeqc in PhreeqcMatlab
clc

phreeqc_rm = PhreeqcSingleCell('methanation_zvi_input.phr' , 'llnl.dat');

h_out = phreeqc_rm.GetSelectedOutputHeadings(1);

dt = 10000; %[s]
t_end = 1000*dt; %[s]
t = dt:dt:t_end;
n_data = length(dt:dt:t_end);
reactants = zeros(n_data, length(h_out));

for i = 1:length(t)
    phreeqc_rm.RM_SetTime(t(i));
    phreeqc_rm.RM_SetTimeStep(dt);

    status = phreeqc_rm.RM_RunCells();
    % c_out = phreeqc_rm.GetConcentrations();

    reactants(i, :) = phreeqc_rm.GetSelectedOutput(1);
    disp(i/length(t)*100)
    if reactants(i,1)<=1e-2
        break;
    end
end
figure(1)
yyaxis left
plot(t, reactants(:,1))
ylabel(h_out{1})
yyaxis right
plot(t, reactants(:,2))
ylabel(h_out{2})

figure(2)
yyaxis left
plot(t, reactants(:, 4));
ylabel(h_out{4})
yyaxis right
plot(t, reactants(:, 7));
ylabel(h_out{7})

figure(3)
plot(t, reactants(:, end));
xlabel('time[s]')
ylabel('zvi mol')

figure(4)
plot(t, reactants(:, 8));
ylabel(h_out{8})

save('methanation_zvi.mat', 'reactants');
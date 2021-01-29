function [new_elec,current_demand,surplus_elec,shortage_elec] = extract_surplus(d, future_coef)
% extract the surplus, deficit, and other relevant values for data d
% future_coef is a multiplier for that is used to calculate the future
% wind electricity production in Denmark
% Example
% date1 = '2021-01-01';
% date2 = '2021-01-20';
% extract_surplus(date1, date2);

d.HourDK = datetime(d.HourDK, 'Format', 'yyyy-MM-dd''T''HH:mm:ss');

new_elec = future_coef*(d.OnshoreWindPower+d.OffshoreWindPower);
surplus_elec = new_elec-d.TotalLoad;
new_elec(new_elec>d.TotalLoad) = d.TotalLoad(new_elec>d.TotalLoad);
surplus_elec(surplus_elec<0) = 0;
shortage_elec = d.TotalLoad-new_elec;
current_demand = d.TotalLoad;

end


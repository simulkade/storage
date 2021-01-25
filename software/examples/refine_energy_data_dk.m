% this scripts refines the electricity data of the last two years 
% for Denmark downloaded from https://www.energidataservice.dk/tso-electricity/electricitybalancenonv

% Data represents the overall balance of consumption, production, import and export of electricity in an area.
% Production is divided into main production types. Data is based on online power measurements for SCADA, and therefore with some delay.
% Gross consumption = sum of production + sum of exchange to connected areas. A positive exchange is import of electricity, while a negative is export.
% The total production is the sum of all production types.

% All the units are MWh/h (power)
% MWh/h = MW



% read the data into a table:
d = readtable('electricitybalancenonv.csv');
d = flip(d); % to have them in increasing date order

% convert the DK hour format from text to datetime
d.HourDK = datetime(d.HourDK, 'Format', 'yyyy-MM-dd''T''HH:mm:ss');

% get rid of some of the columns
d.HourUTC = [];
d.PriceArea = [];

d1 = d(1:2:end-1, :);
d2 = d(2:2:end, :);
d_new = d1; % a temp variable

% create a new table by adding odd and even rows
% data is for two parts of Denmark and we need overal values
f = d.Properties.VariableNames;
for i = 1:length(f)
    if strcmp(f{i}, 'HourDK')
        continue
    end
    d_new{:, f{i}} = d1{:, f{i}}+d2{:, f{i}};
end

% write it to a file
writetable(d_new, 'elec_data_dk.csv');


% plot(d_new.HourDK, d_new.OffshoreWindPower)

% Now extract the data from January 2021
% subset = [datenum(2021, 1,1); datenum(2021, 1, 20)];

% date1 = datenum(2021, 1, 10);
% date2 = datenum(2021, 1, 25);

date1 = datenum(2020, 1, 1);
date2 = datenum(2020, 2, 25);

date_all = datenum(d_new.HourDK);

ind1 = find(date_all==date1, 1);
ind2 = find(date_all==date2, 1);

d_recent = d_new(ind1:ind2, :);

plot(d_recent.HourDK, d_recent.TotalLoad)
hold on
plot(d_recent.HourDK, d_recent.OffshoreWindPower+d_recent.OnshoreWindPower)
hold off

% Notes: 
% 1- write a function to obtain data from a certain period
% 2- create the simple cases for the calculation of energy balance
% 3- work on the subsurface simulator
% 4- finish the equipment simulations and put something together for the report
% approach: there are still lots of things left to do; focus!





function average_val = average_values(t, val)
%average_values returns the average value of val by integrating over time t
d_date = datetime(t, 'Format', 'yyyy-MM-dd''T''HH:mm:ss');
date_num = datenum(d_date);
average_val = trapz(date_num, val)/(date_num(end)-date_num(1));
end


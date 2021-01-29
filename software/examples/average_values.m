function average_val = average_values(t, val)
%average_values returns the average value of val by integrating over time t
d.HourDK = datetime(d.HourDK, 'Format', 'yyyy-MM-dd''T''HH:mm:ss');
date_num = datenum(d.HourDK);
average_val = trapz(date_num, val)/(t(end)-t(1));
end


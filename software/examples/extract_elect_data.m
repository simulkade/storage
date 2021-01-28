function d = extract_elect_data(begin_date, end_date, file_name)
    % d = extract_elec_data(beging_date, end_date)
    % extracts data from the cleaned up data file of electricity production and consumption in denmark
    % data is from the last two years, cleaned up  a bit to give the overal electricity demand
    % and supply of Denmark
    % d is a table of data for the period between begin_date and end_date
    % dates have the following format
    % Example:
    % date1 = '2021-01-15';
    % or
    % date1 = '2021-Jan-15';
    % example use:
    % date1 = '2020-June-13';
    % date2 = '2021-Jan-25';
    % d = extract_elect_data(date1, date2)

    if nargin<3
        file_name = 'elec_data_dk.csv';
    end
    
    d_base = readtable(file_name);

    date1 = datenum(begin_date);
    date2 = datenum(end_date);

    date_all = datenum(datetime(d_base.HourDK, 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
    % date_all = datenum(d_base.HourDK);

    ind1 = find(date1<=date_all, 1)
    ind2 = find(date2<=date_all, 1)
    
    d = d_base(ind1:ind2, :);

end

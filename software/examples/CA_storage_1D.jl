using DataFrames, CSV, Dates

file_name = "elec_data_dk.csv"

# load the data file
df = DataFrame(CSV.File(file_name))

# separate the relevant section
date_begin = DateTime("2020-10-01")
date_end = DateTime("2021-01-20")
d_long = df[date_begin .<= df.HourDK .<= date_end, :]
d_long.HourDK .+= Dates.Year(30)

# extract the surplus and demand in the future
future_coef = 4.0
demand_coef = 1.5
new_elec = future_coef*(d_long.OnshoreWindPower+d_long.OffshoreWindPower)
current_demand = demand_coef*d_long.TotalLoad
surplus_elec = new_elec-current_demand
new_elec[new_elec .> current_demand] .= current_demand[new_elec .> current_demand]
surplus_elec[surplus_elec .< 0] .= 0.0
shortage_elec = current_demand-new_elec

# set the end of storage
end_storage = "2020-11-30"
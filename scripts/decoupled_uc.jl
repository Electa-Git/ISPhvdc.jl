# Define all packages that are needed
using Pkg
using ISPhvdc
using PowerModels
using PowerModelsACDC
using CbaOPF
# using NEM2000synthetic
using JuMP
using Ipopt
using Plots
using PlotlyJS
using PowerPlots
using DataFrames
using CSV
using Gurobi
using StatsPlots
using StatsBase
using Statistics
using JSON
using Plots
using HiGHS

# Create short hands for the most important ones
const _PM = PowerModels
const _PMACDC = PowerModelsACDC
# const _SNEM = NEM2000synthetic
const _PP = PowerPlots
const _SB = StatsBase
const _ISP = ISPhvdc
cd("/Users/hergun/.julia/dev/ISPhvdc")
###########################################################
################## INPUT SECTION ##########################
# SELECT SCENARIO
# scenario ∈ {"2022 ISP Hydrogen Superpower", "2022 ISP Progressive Change", "2022 ISP Slow Change", "2022 ISP Step Change"}
scenario = "2022 ISP Step Change"
# select climate year ∈ [2024:2051]
year = 2034
# Number generator contingencies per area 1 = NSW, 2 = VIC, 3 = QLD, 4 = SA, 5 = TAS.
generator_contingencies = [50 30 50 30 20] 
# You can choose select certain hours or a full year for the analysis: 
# selected_hours = Dict{String, Any}("hour_range" => start hour:end hour)
# selected_hours = Dict{String, Any}("all")
selected_hours = Dict{String, Any}("hour_range" => 1:2:48) #8737:2:8784
# State if data ISP should be downloaded, only necessary for the first time, takes about 3 minutes!
download_data = false
# State if circiuts and parallel lines should be merged:
merge_parallel_lines = true
# Assign solvers
dc_solver =  JuMP.optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 0.5e-2, "TimeLimit" => 900) # "mip_focus" => 3 https://www.gurobi.com/documentation/current/refman/method.html#parameter:Method 
dn_res_factor = 0.0
t_fcrd = 6.0
t_fcr = 1.0
t_hvdc = 0.5
extension = "_rd_reserves_fdb"
############ END INPUT SECTION ##############################
#############################################################


#############################################################
#################### START METHODOLOGY ######################
if download_data == true
    _ISP.download_isp_data()
end

# Optimisation settings for CbaOPF.jl
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "uc_reserves" => true)

# Test case data
data_folder = joinpath("data")
data_file_hvdc = "nem_2300bus_hvdc.m" #"nem_2300bus_thermal_limits_gen_costs_hvdc_v1.m"

# Get grid data from the NEM 2000 bus model m-file 
data = _PM.parse_file(data_folder*"/"*data_file_hvdc)
# Process data to fit into PMACDC model
_PMACDC.process_additional_data!(data)
# Delete DC lines which have been modelled as AC lines
_ISP.fix_hvdc_data_issues!(data)
# Assign buses to states
_ISP.add_area_dict!(data)
# Extend data model with flexible demand, to be able to do demand shedding if required
_ISP.add_demand_data!(data)
# Aggregate demand data per state to modulate with hourly traces
_ISP.aggregate_demand_data!(data)
# Get demand traces for selected year, for each state
total_demand_series = _ISP.get_demand_data(scenario, year)
average_demand_per_state = Dict{String, Any}([state => mean(timeseries) for (state, timeseries) in total_demand_series])
# Create demand series for whole Australia for plotting and inspection
aggregated_demand_series = sum([ts for (state, ts) in total_demand_series])
# Plot the demand series
Plots.plot(aggregated_demand_series[1:48])
# Get the load substractors, e.g., distirubted generation data per state
dn_demand_series  = _ISP.get_dn_demand_data(scenario, year)
# Get installed generator info, e.g. installed generation type and capacity from the ISP data
generator_info = _ISP.get_generator_information()
# Get RES time series, e.g. traces from the ISP data
pv, wind = _ISP.get_res_timeseries(year)
# Aggregate timeseries to obtain one profile for existing RES
pv_series, count_pv = _ISP.aggregate_res_timeseries(pv, generator_info, "Solar")
wind_series, count_wind = _ISP.aggregate_res_timeseries(wind, generator_info, "Wind")
# Aggregate timeseries to obtain one profile for renewable energy zones (REZ)
pv_rez = _ISP.make_rez_time_series(pv)
wind_rez = _ISP.make_rez_time_series(wind)

# Get generation capacity of REZ and the grid extensions and update grid data
rez_capacities = _ISP.get_rez_capacity_data(scenario, year)
rez_connections = _ISP.get_rez_grid_extensions()
_ISP.add_rez_and_connections!(data, rez_connections, rez_capacities)

# fix data issues, e.g. putting generation cost in € / pu:
_ISP.fix_data!(data)

if merge_parallel_lines == true
    _ISP.merge_parallel_lines(data)
end

# Select hours
hours = _ISP.select_hours(year, selection = selected_hours)

fmin = 49.5
data["frequency_parameters"] = Dict{String, Any}()
data["frequency_parameters"]["fmin"] = fmin
data["frequency_parameters"]["f0"] = 50.0
data["frequency_parameters"]["fdb"] = 0.1
data["frequency_parameters"]["fmax"] =  data["frequency_parameters"]["f0"] + ((data["frequency_parameters"]["f0"] - fmin))
data["frequency_parameters"]["t_fcr"] = t_fcr
data["frequency_parameters"]["t_fcrd"] = t_fcrd
data["frequency_parameters"]["uc_time_interval"] = 1.0 # hours
data["frequency_parameters"]["delta_fss"] = 0.2


h = join([hours[1],"_", hours[end]])
y = "$year"
fn = joinpath("results", scenario, y, h, join(["branch_expansion_",dn_res_factor,"_",hours[end],".json"]))
expansion = Dict{String, Any}()
open(fn) do f
    dicttxt = read(f,String)  # file information to string
    global expansion = JSON.parse(dicttxt)  # parse and transform data
end

for (b, branch) in data["branch"]
    if haskey(expansion["solution"]["nw"]["1"]["branch"], b)
        branch["rate_a"] = branch["rate_a"] + (expansion["solution"]["nw"]["1"]["branch"][b]["delta_cap"] * 1.0)
    else
        branch["rate_a"] = branch["rate_a"] * 1.0
    end
    branch["angmin"] = - (branch["rate_a"] / max(1e-12, branch["br_x"]))
    branch["angmax"] =  (branch["rate_a"] / max(1e-12, branch["br_x"]))
end

# Function to write generator & converter info:
_ISP.generator_uc_data!(data, fcr_cost = 20, droop_fac = 1)
_ISP.converter_uc_data!(data, t_hvdc = t_hvdc, ffr_cost = 20)
_ISP.define_tie_lines!(data)

fmin = 49.0:0.1:49.7
droop = true

data_dict = Dict()
data_dict["data"] = data
data_dict["total_demand_series"] = total_demand_series
data_dict["dn_demand_series"] = dn_demand_series
data_dict["pv_series"] =  pv_series
data_dict["wind_series"] = wind_series
data_dict["wind_rez"] = wind_rez
data_dict["pv_rez"] = pv_rez
data_dict["generator_contingencies"] = generator_contingencies
data_dict["no_dc_cont"] = false
data_dict["dn_res_factor"] = dn_res_factor
data_dict["p2p"] = false


@time mn_data = _ISP.multi_network_uc_data_no_cont(data, total_demand_series, dn_demand_series, pv_series, wind_series, pv_rez, wind_rez, hours; dn_res_factor = dn_res_factor)

@elapsed r_uc = CbaOPF.solve_uc(mn_data, _PM.DCPPowerModel, dc_solver, setting = s, multinetwork = true)

objective_dc, objective_no_dc, time_dc, time_no_dc, time_opt_dc, time_opt_no_dc = _ISP.batch_fsrdopf(data_dict, fmin, dc_solver, scenario, year, hours, h, r_uc; droop = true, extension = extension)

_ISP.get_and_plot_objective_value(fmin, scenario, "$year", h; extension = extension)

year_ = "$year"

_ISP.plot_largest_continegncy(scenario, year_, h, "49.0", data_dict, 1:length(hours); extension = extension)
_ISP.plot_largest_continegncy(scenario, year_, h, "49.5", data_dict, 1:length(hours); extension = extension)

_ISP.plot_calculation_time(fmin, 49.0, scenario, year_, h; extension = extension)
_ISP.plot_calculation_time(fmin, 49.5, scenario, year_, h; extension = extension)
_ISP.plot_hvdc_contribution(data_dict, 49.0, scenario, year_, h, 1:length(hours); extension = extension)
_ISP.plot_hvdc_contribution(data_dict, 49.5, scenario, year_, h, 1:length(hours); extension = extension)


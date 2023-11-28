# Define all packages that are needed
using Pkg
using ISPhvdc
using PowerModels
using PowerModelsACDC
using CbaOPF
using NEM2000synthetic
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

# Create short hands for the most important ones
const _PM = PowerModels
const _PMACDC = PowerModelsACDC
const _SNEM = NEM2000synthetic
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
# You can choose select certain hours or a full year for the analysis: 
# selected_hours = Dict{String, Any}("hour_range" => start hour:end hour)
# selected_hours = Dict{String, Any}("all")
selected_hours = Dict{String, Any}("hour_range" => 368:1:368)
# State if data ISP should be downloaded, only necessary for the first time, takes about 3 minutes!
download_data = false
# Select OPF method opf ∈ {"AC", "DC", "LPAC", "SOC"}
opf = "DC"
# State if circiuts and parallel lines should be merged:
merge_parallel_lines = true
# Assign solvers
ac_solver =  JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000, "print_level" => 0, "hsllib" => "/Users/hergun/IpoptMA/lib/libhsl.dylib", "linear_solver" => "ma27")
dc_solver =  JuMP.optimizer_with_attributes(Gurobi.Optimizer) #  https://www.gurobi.com/documentation/current/refman/method.html#parameter:Method 
lpac_solver =  JuMP.optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0)
soc_solver =  JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000, "print_level" => 0, "hsllib" => "/Users/hergun/IpoptMA/lib/libhsl.dylib", "linear_solver" => "ma27")

dn_res_factor = 0.0
############ END INPUT SECTION ##############################
#############################################################


#############################################################
#################### START METHODOLOGY ######################
if download_data == true
    _ISP.download_isp_data()
end

# Optimisation settings for CbaOPF.jl
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true)

# Test case data
data_folder = joinpath("data")
data_file_hvdc = "nem_2300bus_hvdc.m"

# Get grid data from the NEM 2000 bus model m-file 
data_hvdc = _PM.parse_file(data_folder*"/"*data_file_hvdc)
# Process data to fit into PMACDC model
_PMACDC.process_additional_data!(data_hvdc)
# Delete DC lines which have been modelled as AC lines
_ISP.fix_hvdc_data_issues!(data_hvdc)
# Assign buses to states
_ISP.add_area_dict!(data_hvdc)
# Extend data model with flexible demand, to be able to do demand shedding if required
_ISP.add_demand_data!(data_hvdc)
# Aggregate demand data per state to modulate with hourly traces
_ISP.aggregate_demand_data!(data_hvdc)

# Get demand traces for selected year, for each state
total_demand_series = _ISP.get_demand_data(scenario, year)
average_demand_per_state = Dict{String, Any}([state => mean(timeseries) for (state, timeseries) in total_demand_series])
# Create demand series for whole Australia for plotting and inspection
aggregated_demand_series = sum([ts for (state, ts) in total_demand_series])
# Plot the demand series
Plots.plot(aggregated_demand_series')
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
_ISP.add_rez_and_connections!(data_hvdc, rez_connections, rez_capacities)

# Make copy of grid data for hourly calculations
opf_data = deepcopy(data_hvdc)


# fix data issues, e.g. putting generation cost in € / pu:
_ISP.fix_data!(opf_data)

if merge_parallel_lines == true
    _ISP.merge_parallel_lines(opf_data)
end

# Create empty arrays
total_demand = []
total_gen_capacity = []
# Make copy of grid data
hourly_data = deepcopy(opf_data)

for (b, branch) in hourly_data["branch"]
    branch["delta_cap_max"] = branch["rate_a"] * 2 # for testing.....
    branch["capacity_cost"] = 150e4 * hourly_data["baseMVA"] / (25 * 8760) # for testing, update with more realistic numbers.....
end

for (l, load) in hourly_data["load"]
    load["cost_curt"] = 1e6
end

fmin = 49.5
hourly_data["frequency_parameters"] = Dict{String, Any}()
hourly_data["frequency_parameters"]["fmin"] = fmin
hourly_data["frequency_parameters"]["f0"] = 50.0
hourly_data["frequency_parameters"]["fmax"] =  hourly_data["frequency_parameters"]["f0"] + ((hourly_data["frequency_parameters"]["f0"] - fmin))
hourly_data["frequency_parameters"]["t_fcr"] = 0.1
hourly_data["frequency_parameters"]["t_fcrd"] = 6.0
hourly_data["frequency_parameters"]["uc_time_interval"] = 1.0 # hours

# Select hours
hours = _ISP.select_hours(year, selection = selected_hours)

_ISP.generator_uc_data!(hourly_data)

mn_data =  _ISP.multi_network_capacity_data(hourly_data, total_demand_series, dn_demand_series, pv_series, wind_series, pv_rez, wind_rez, hours, dn_res_factor = dn_res_factor )

result = CbaOPF.solve_nodal_tnep(mn_data, _PM.DCPPowerModel, dc_solver; multinetwork = true, setting = s)

for (b, branch) in result["solution"]["nw"]["1"]["branch"]
    if branch["delta_cap"] .!== 0.0
        print(b, " ", branch["delta_cap"], "\n")
    end
end

ls = sum( [sum([load["pcurt"] for (l, load) in result["solution"]["nw"][nw]["load"]]) for (nw, network) in result["solution"]["nw"]])

print("load shedding = " , ls*100, "\n")


y = "$year"
h = join([hours[1],"_", hours[end]])
if !isdir(joinpath("results", scenario))
    mkdir(joinpath("results", scenario))
end
if !isdir(joinpath("results", scenario, y))
    mkdir(joinpath("results", scenario, y))
end
if !isdir(joinpath("results", scenario, y, h))
    mkdir(joinpath("results", scenario, y, h))
end

filename = joinpath("results", scenario, y, h, join(["branch_expansion_",dn_res_factor,"_",hours[end],".json"]))
json_string = JSON.json(result)
open(filename,"w") do f
write(f, json_string)
end
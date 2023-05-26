# Define all packages that are needed
using Pkg
using ISPhvdc
using PowerModels
using PowerModelsACDC
using CbaOPF
using NEM_2300bus
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

# Create short hands for the most important ones
const _PM = PowerModels
const _PMACDC = PowerModelsACDC
const _SNEM = NEM_2300bus
const _PP = PowerPlots
const _SB = StatsBase
const _ISP = ISPhvdc
###########################################################
################## INPUT SECTION ##########################

# SELECT SCENARIO
# scenario ∈ {"2022 ISP Hydrogen Superpower", "2022 ISP Progressive Change", "2022 ISP Slow Change", "2022 ISP Step Change"}
scenario = "2022 ISP Step Change"
# select climate year ∈ [2024:2051]
year = 2024
# You can choose select certain hours or a full year for the analysis: 
# selected_hours = Dict{String, Any}("hour_range" => start hour:end hour)
# selected_hours = Dict{String, Any}("all")
selected_hours = Dict{String, Any}("hour_range" => 1:12)
# State if data ISP should be downloaded, only necessary for the first time, takes about 3 minutes!
download_data = false
# Select OPF method opf ∈ {"AC", "DC", "LPAC"}
opf = "LPAC"
# Assign solvers
ac_solver =  JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000, "print_level" => 0, "hsllib" => "/Users/hergun/IpoptMA/lib/libhsl.dylib", "linear_solver" => "ma27")
dc_solver =  JuMP.optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0)


############ END INPUT SECTION ##############################
#############################################################


#############################################################
#################### START METHODOLOGY ######################
if download_data == true
    _ISP.download_isp_data()
end

# Optimisation settings for CbaOPF.jl
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "fix_cross_border_flows" => false)

# Test case data
data_folder = joinpath(dirname(pathof(NEM_2300bus)),"..","data/matpower")
data_file_hvdc = "nem_2300bus_thermal_limits_gen_costs_hvdc.m"

# Get grid data from the NEM 2000 bus model m-file 
data_hvdc = _PM.parse_file(data_folder*"/"*data_file_hvdc)
# Process data to fit into PMACDC model
_PMACDC.process_additional_data!(data_hvdc)
# Delete DC lines which have been modelled as AC lines
_SNEM.fix_hvdc_data_issues!(data_hvdc)
# Assign buses to states
_SNEM.add_area_dict!(data_hvdc)
# Extend data model with flexible demand, to be able to do demand shedding if required
_ISP.add_demand_data!(data_hvdc)

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
# Aggregate demand data per state to modulate with hourly traces
_ISP.aggregate_demand_data!(opf_data)

# Create empty arrays
total_demand = []
total_gen_capacity = []
# Make copy of grid data
hourly_data = deepcopy(opf_data)
# Select hours
hours = _ISP.select_hours(year, selection = selected_hours)

# Optionally, determine hosting capacity for all nodes
# hosting_capacity = calculate_hosting_capacity(opf_data, total_demand_series, dn_demand_series, hours)

# Create dictionaries for inspection of results
pf = Dict{String, Any}(["$hour" => zeros(length(opf_data["branch"])) for hour in hours])
pfdc = Dict{String, Any}(["$hour" => zeros(length(opf_data["branchdc"])) for hour in hours])
pcurt = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])
pd = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])
pgmax = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])
pg = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])

# Run hourly OPF calcuations
for hour in hours
    # Write hourly pv, wind and demand traces into opf data
    _ISP.prepare_hourly_opf_data!(hourly_data, opf_data, total_demand_series, average_demand_per_state, pv_series, wind_series, pv_rez, wind_rez, hour)
    # calculate some charactersitics for result inspecttion
    pd_max = sum([load["pd"] for (l, load) in opf_data["load"]]) * hourly_data["baseMVA"]
    pg_max = sum([gen["pmax"] for (g, gen) in opf_data["gen"]]) * hourly_data["baseMVA"]
    pdh_max = sum([load["pd"] for (l, load) in hourly_data["load"]]) * hourly_data["baseMVA"]
    pgh_max = sum([gen["pmax"] for (g, gen) in hourly_data["gen"]]) * hourly_data["baseMVA"]
    # print charactersiticss
    print("Grid Data: Total demand = ", pd_max, " MW, Total generation = ", pg_max, " MW","\n")
    print("Hour ", hour,": Total demand = ", pdh_max, " MW, Total generation = ", pgh_max, " MW","\n")

    push!(total_demand, sum([load["pd"] for (l, load) in hourly_data["load"]]))
    push!(total_gen_capacity, sum([gen["pmax"] for (g, gen) in hourly_data["gen"]]))
    # Solve OPF
    if opf == "AC"
        opf_result = CbaOPF.solve_cbaopf(hourly_data, _PM.ACPPowerModel, ac_solver, setting = s)
    elseif opf == "LPAC"
        opf_result = CbaOPF.solve_cbaopf(hourly_data, _PM.LPACCPowerModel, dc_solver, setting = s)
    elseif opf == "DC"
        opf_result = CbaOPF.solve_cbaopf(hourly_data, _PM.DCPPowerModel, dc_solver, setting = s)
    end
    # Write out general information
    print("Hour: ", hour, " -> ", opf_result["termination_status"], " in ", opf_result["solve_time"], " seconds.", "\n")
    # calculate and print some more information base on results
    if haskey(opf_result["solution"], "load")
        pd["$hour"] = [hourly_data["load"]["$l"]["pd"] for l in sort(parse.(Int, collect(keys(opf_result["solution"]["load"]))))]
        pgmax["$hour"] = [hourly_data["gen"]["$g"]["pmax"] for g in sort(parse.(Int, collect(keys(opf_result["solution"]["gen"]))))]
        pcurt["$hour"] = [opf_result["solution"]["load"]["$l"]["pcurt"] for l in sort(parse.(Int, collect(keys(opf_result["solution"]["load"]))))]
        pcurt_tot = sum([opf_result["solution"]["load"]["$l"]["pcurt"] for l in sort(parse.(Int, collect(keys(opf_result["solution"]["load"]))))]) * hourly_data["baseMVA"]
        pf["$hour"] = [opf_result["solution"]["branch"]["$b"]["pf"] for b in sort(collect(parse.(Int, keys(opf_result["solution"]["branch"]))))] ./ [hourly_data["branch"]["$b"]["rate_a"] for b in sort(parse.(Int, collect(keys(opf_result["solution"]["branch"]))))]
        pfdc["$hour"] = [opf_result["solution"]["branchdc"]["$b"]["pf"] for b in sort(parse.(Int, collect(keys(opf_result["solution"]["branchdc"]))))] ./ [hourly_data["branchdc"]["$b"]["rateA"] for b in sort(parse.(Int, collect(keys(opf_result["solution"]["branchdc"]))))]
        pg["$hour"] = [opf_result["solution"]["gen"]["$g"]["pg"] for g in sort(parse.(Int, collect(keys(opf_result["solution"]["gen"]))))]
    else
        pcurt_tot = 0
    end
    print("Total curtailed load = ", pcurt_tot, " MW, ", pcurt_tot / pdh_max * 100,"%", "\n")
end


##### HERE there is some loose code for inspecting results


# bus_loads = Dict{String, Any}([b => [] for (b, bus) in opf_data["bus"]])
# for (l, load) in opf_data["load"]
#     lb = load["load_bus"]
#     push!(bus_loads["$lb"], parse(Int, l))
# end

# pcurt_avg = sum([pcurt[h] for (h, hour) in pcurt]) ./ length(pcurt)
# pd_avg = sum([pd[h] for (h, hour) in pd]) ./ length(pd)
# Plots.plot(pcurt_avg .* opf_data["baseMVA"])

# pmax_avg = sum([pgmax[h] for (h, hour) in pgmax]) ./ length(pgmax)
# pg_avg = sum([pg[h] for (h, hour) in pg]) ./ length(pg)
# Plots.scatter(pmax_avg .- pg_avg)


# pf_avg = sum([pf[h] for (h, hour) in pf]) ./ length(pf)
# Plots.plot(pf_avg) 

# for idx in 1:length(pf_avg)
#     if pf_avg[idx] > 0.8
#         fbus = opf_data["branch"]["$idx"]["f_bus"]
#         tbus = opf_data["branch"]["$idx"]["t_bus"]
#         loadsf = bus_loads["$fbus"]
#         loadst = bus_loads["$tbus"]
#         lsf = sum(pcurt_avg[loadsf]) * opf_data["baseMVA"]
#         lst = sum(pcurt_avg[loadst]) * opf_data["baseMVA"]
#         print("Branch ", idx, " fbus = ", fbus, ", tbus = ", tbus, " -> load shedding fbus = ", lsf," MW, load shedding tbus = ", lst, " MW", "\n" )
#     end 
# end

# for idx in 1:length(pcurt_avg)
#     if pcurt_avg[idx] > 0.3
#         print("Load shedding in load ", idx, " connected at bus ", opf_data["load"]["$idx"]["load_bus"], " = ", pcurt_avg[idx] * opf_data["baseMVA"], " MW","\n")
#     end
# end


# for idx in 1:length(pmax_avg)
#     if (pmax_avg[idx] - pg_avg[idx]) > 5
#         id = gen_ids[idx]
#         print("Generation curtailment in generator ", id, " connected at bus ", opf_data["gen"]["$id"]["gen_bus"], " = ", (pmax_avg[idx] - pg_avg[idx]) * opf_data["baseMVA"], " MW","\n")
#     end
# end

# Plots.plot(total_demand .* data_hvdc["baseMVA"], label = "Total demand in MW")
# Plots.plot!(total_gen_capacity .* data_hvdc["baseMVA"], label = "Total available generation in MW")
# Plots.plot((total_demand .- total_gen_capacity).* data_hvdc["baseMVA"] , label = "ΔP in MW")

# for (b, branch) in opf_data["branch"]
#     if branch["f_bus"] == 2013 || branch["t_bus"] == 2013 
#         print(b, " ", branch["rate_a"], " fbus = ", branch["f_bus"], " tbus = ", branch["t_bus"], "\n")
#     end
# end
# StatsPlots.boxplot!(hosting_capacity["1402"]')

# a = Statistics.quantile(hosting_capacity["1402"], 0.95)

# for (l, load) in opf_result["1"]["1"]["36"]["solution"]["load"]
#     if load["pcurt"] >= 1e-5
#      print("Load ", l, " -> curtailed demand = ", load["pcurt"]*100 , " MW", "\n")
#     end
# end 


# result = opf_result["1"]["1"]["2"]
# time_stamp = 2
# @time hourly_data = prepare_hourly_opf_data!(hourly_data, opf_data, total_demand_series, dn_demand_series, pv_series, wind_series, time_stamp)
# load_orig = [opf_data["load"]["$l"]["pd"] for l in sort(collect(keys(opf_data["load"])))]
# load_hour = [hourly_data["load"]["$l"]["pd"] for l in sort(collect(keys(hourly_data["load"])))]
# load_flex = [result["solution"]["load"]["$l"]["pflex"] for l in sort(collect(keys(hourly_data["load"])))]
# Plots.plot(load_hour, load_flex, seriestype=:scatter, markersize = 5)
# Plots.plot!(-1:10, -1:10)

# print("Total curtailed load = ", sum([result["solution"]["load"]["$l"]["pcurt"] for l in sort(collect(keys(hourly_data["load"])))]) * data_hvdc["baseMVA"], "\n")

# pg_orig = [opf_result_grid["solution"]["gen"]["$g"]["pg"] for g in sort(collect(keys(opf_result_grid["solution"]["gen"])))]
# pg_res = [result["solution"]["gen"]["$g"]["pg"] for g in sort(collect(keys(result["solution"]["gen"])))]
# Plots.plot(pg_orig, pg_res, seriestype=:scatter, markersize = 2)
# Plots.plot!(1:10, 1:10)

# result = opf_result["1"]["1"]["2"]
# pf_orig = [opf_result_grid["solution"]["branch"]["$b"]["pf"]/data_hvdc["branch"]["$b"]["rate_a"] for b in sort(collect(keys(opf_result_grid["solution"]["branch"])))]
# pf_res = [result["solution"]["branch"]["$b"]["pf"]/data_hvdc["branch"]["$b"]["rate_a"] for b in sort(collect(keys(result["solution"]["branch"])))]
# Plots.plot(pf_orig, pf_res, seriestype=:scatter, markersize = 2)
# Plots.plot!(0:1, 0:1)
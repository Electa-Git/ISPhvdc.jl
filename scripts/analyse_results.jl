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
######################################
scenario = "2022 ISP Step Change"
year = "2034"
hours = "1_47"
#hours = "8737_8783"
fmin = 49.0:0.1:49.8

input_data = _ISP.load_input_data(scenario, year, hours)
_ISP.plot_system_information(input_data, scenario, year, hours)

# print("Savings = ", sum([objective for (o, objective) in objective_no_dc]) - sum([objective for (o, objective) in objective_dc]), "\n")

objective_dc, objective_no_dc  = _ISP.get_and_plot_objective_value(fmin, scenario, year, hours)
print([ob - objective_dc[o] for (o, ob) in objective_no_dc])
_ISP.plot_calculation_time(fmin, scenario, year, hours)


fmin_ = 49.7
_ISP.plot_load_shedding(input_data, fmin_, scenario, year, hours)
_ISP.plot_total_inertia(input_data, fmin_, scenario, year, hours)
_ISP.plot_tie_line_flows(input_data, fmin_, scenario, year, hours)
# _ISP.plot_dc_flows(input_data, fmin_, scenario, year, hours)
_ISP.plot_res_generation_and_curtailment(input_data, fmin_, scenario, year, hours)
_ISP.plot_hvdc_contribution(input_data, fmin, scenario, year, hours)


fn = joinpath("results",scenario, year, hours, join(["f",49.8,"_with_dc.json"]))
result_dc = Dict{String, Any}()
open(fn) do f
dicttxt = read(f,String)  # file information to string
    global result_dc = JSON.parse(dicttxt)  # parse and transform data
end


# # for (l, load) in result_dc["solution"]["nw"]["61"]["load"]
# #     if load["pcurt"] >= 1e-9
# #     print(l, ": ", load["pcurt"], "\n")
# #     end
# # end

# for (c, conv) in result_dc["solution"]["nw"]["1"]["convdc"]
#     print(c, " ", conv["pdc"], "\n")
# end
# # for n in sort(parse.(Int, keys(input_data["nw"])))
#     nw = input_data["nw"]["$n"]
#     number_of_contingencies = input_data["number_of_contingencies"]
#     if mod(n, number_of_contingencies) == 0
#         hour_id = Int(n - number_of_contingencies + 1)
#     else
#         hour_id = Int(n - mod(n, number_of_contingencies) + 1)
#     end
#     if !isnothing(nw["contingency"]["gen_id"])
#         gen_id = nw["contingency"]["gen_id"]
#         print(hour_id, " ", gen_id, ": ", nw["gen"]["$gen_id"]["name"], " ", nw["gen"]["$gen_id"]["type"], ", pmax = ",input_data["nw"]["$hour_id"]["gen"]["$gen_id"]["pmax"],  "\n")
#     end
# end

# for i in a
#     print(i, " ", input_data["nw"]["61"]["gen"]["$i"]["name"], " ", input_data["nw"]["61"]["gen"]["$i"]["type"], " ",input_data["nw"]["61"]["gen"]["$i"]["pmax"], "\n")
# end

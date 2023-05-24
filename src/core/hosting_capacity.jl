# This function calculated the theoretic hosting capacity of each node defined as:
# HC = p(∑(all loads connected to the bus)_t + ∑(all line ratings connected to the bus)_t)
# Inputs:
# - grid_data: PowerModelsACDC dictionary
# - demand_series: half hourly demand series per state
# - dn_demand_series: half hourly load subtractor values indicating distributed generation
# - hours: selected hours to run the analysis
# - percentile ∈ [0, 1], e.g. p() = 0.5, 0.9, 0.95, ....
# - hourly: This input gives the hosting capacity "time series" as output instead of its percentile value

function calculate_hosting_capacity(grid_data, demand_series, dn_demand_series, hours; percentile = 0.5, hourly = true)
    hosting_capacity = Dict{String, Any}([b => zeros(1, length(hours)) for (b, bus) in grid_data["bus"]])

    for (b, bus) in grid_data["bus"]
        print("Processing bus: ", b, "\n")
        calculate_demand!(grid_data, demand_series, dn_demand_series, hours, b, hosting_capacity)
        line_capacities = calculate_line_capacities(grid_data, bus)
        hosting_capacity[b] = hosting_capacity[b] .+ line_capacities
    end


    return hosting_capacity
end


function calculate_line_capacities(grid_data, bus)
    line_cap = 0

    for (b, branch) in grid_data["branch"]
        if branch["f_bus"] == bus["index"] || branch["t_bus"] == bus["index"]
            line_cap = line_cap + branch["rate_a"]
        end
    end

    for (c, conv) in grid_data["convdc"]
        if conv["busac_i"] == bus["index"]
            line_cap = line_cap + conv["Pacrated"]
        end
    end

    return line_cap
end

function calculate_demand!(grid_data, demand_series, dn_demand_series, hours, bus, hosting_capacity)
    for (l, load) in grid_data["load"]
        if load["load_bus"] == grid_data["bus"][bus]["index"]
            h_idx = 1
            for hour in hours
                demand = 0
                area_code = grid_data["bus"][bus]["area"]
                area = grid_data["areas"]["$area_code"]
                area_demand_grid_data = grid_data["aggregated_data"][area]["demand"]
                demand_trace = demand_series[area][hour]
                dn_wind_trace = dn_demand_series["Wind"][area][hour]
                if haskey(dn_demand_series["PV"], area)
                    dn_pv_trace = dn_demand_series["PV"][area][hour]
                else
                    dn_pv_trace = 0
                end

                demand_ratio = (demand_trace[1] + dn_wind_trace[1] + dn_pv_trace[1]) / area_demand_grid_data
                demand = demand + (load["pd"] * demand_ratio)
                hosting_capacity[bus][h_idx] = demand
                h_idx = h_idx + 1
            end
        end
    end
end

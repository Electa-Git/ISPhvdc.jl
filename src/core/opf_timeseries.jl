function prepare_hourly_opf_data!(hourly_data, grid_data, total_demand, dn_demand, pv_series, wind_series, rez_pv, rez_wind, time_stamp)

    hour = time_stamp

    hourly_data["pst"] = Dict{String, Any}()

    for (l, load) in hourly_data["load"]
        load_bus = load["load_bus"]
        area_code = hourly_data["bus"]["$load_bus"]["area"]
        area = hourly_data["areas"]["$area_code"]
        area_demand_grid_data = grid_data["aggregated_data"][area]["demand"]
        demand_trace = total_demand[area][hour]
        demand_ratio = (demand_trace[1] ) / area_demand_grid_data

        # dn_wind_trace = dn_demand["Wind"][area][hour]
        # if haskey(dn_demand["PV"], area)
        #     dn_pv_trace = dn_demand["PV"][area][hour]
        # else
        #     dn_pv_trace = 0
        # end
        #- dn_wind_trace[1] - dn_pv_trace[1]

        if load["pd"] >= 0 
            load["pd"] = grid_data["load"][l]["pd"] * demand_ratio
            load["qd"] = grid_data["load"][l]["qd"] * demand_ratio
        end
    end

    for (g, gen) in hourly_data["gen"]
        trace = 1
        if any(gen["name"] .== keys(rez_pv)) || any(gen["name"] .== keys(rez_wind))
            if any(gen["name"] .== keys(rez_pv)) && gen["type"] == "Solar"
                rez_name = gen["name"]
                trace = rez_pv[rez_name][hour]
            end
            if any(gen["name"] .== keys(rez_wind))  && gen["type"] == "Wind"
                rez_name = gen["name"]
                trace = rez_wind[rez_name][hour]
            end
        elseif gen["gen_status"] == 1
            gen_bus = gen["gen_bus"]
            area_code = hourly_data["bus"]["$gen_bus"]["area"]
            area = hourly_data["areas"]["$area_code"]
            if gen["type"] == "Wind"
                trace = wind_series[area][hour]
            elseif gen["type"] == "Solar"
                if !isempty(pv_series[area])
                    trace = pv_series[area][hour]
                end
            elseif gen["type"] == "VAr support"
                trace = 0
            end
        end
        gen["pmax"] = grid_data["gen"][g]["pmax"] * trace
    end
    if haskey(hourly_data, "aggregated_data")
        delete!(hourly_data, "aggregated_data")
    end
    
    return hourly_data
end

function fix_data!(data)

    # Find isolated buses and put their demand zero
    bus_arcs = Dict{String, Any}([b => [] for (b, bus) in data["bus"]])
    for (b, branch) in opf_data["branch"]
        fbus = branch["f_bus"]
        tbus = branch["t_bus"]
        if haskey(bus_arcs, "$fbus")
            push!(bus_arcs["$fbus"], parse(Int, b))
        end
        if haskey(bus_arcs, "$tbus")
            push!(bus_arcs["$tbus"], parse(Int, b))
        end
    end
    for (c, conv) in opf_data["convdc"]
        cbus = conv["busac_i"]
        push!(bus_arcs["$cbus"], parse(Int, c))
    end

    print(bus_arcs["2013"])

    for (l, load) in data["load"]
        load_bus = load["load_bus"]
        if isempty(bus_arcs["$load_bus"])
            load["pd"] = 0.0
            load["qd"] = 0.0
        end
    end

    return data
end

function aggregate_demand_data!(grid_data)

    grid_data["aggregated_data"] = Dict{String, Any}([area => Dict{String, Any}("demand" => 0) for (key, area) in grid_data["areas"]])

    for (l, load) in grid_data["load"]
        load_bus = load["load_bus"]
        area_code = grid_data["bus"]["$load_bus"]["area"]
        area = grid_data["areas"]["$area_code"]
        grid_data["aggregated_data"][area]["demand"] = grid_data["aggregated_data"][area]["demand"] + load["pd"] * grid_data["baseMVA"]
    end

    return grid_data
end

function get_df_value(data, area, time_stamp; res = false, pv = false)
    if pv == false
        hour = lpad(time_stamp["hour"], 2, "0")
    else
        hour = lpad(time_stamp["hour"], 1, "0")
    end
    if res == false
        month_data =  data[area][data[area][!, "Month"] .== time_stamp["month"], :]
    else
        month_data =  data[area]["values"][data[area]["values"][!, "Month"] .== time_stamp["month"], :]
    end
    trace =  month_data[month_data[!, "Day"] .== time_stamp["day"], hour]

    return trace
end

function find_hours(data, time_stamps)
    df = collect(values(data))[1]


    if !any(data[area][!, "Month"] .== time_stamp["month"])
        output = false
    else
        month_data =  data[area][data[area][!, "Month"] .== time_stamp["month"], :]
        if !any(month_data[!, "Day"] .== time_stamp["day"])
            output = false
        else
            if any(time_stamp["hour"] .> 48)
                output = false
            end
        end
    end

    return output
end

function add_demand_data!(data)
    for (l, load) in data["load"]
        # Superior bound on voluntary load reduction (not consumed power) as a fraction of the total reference demand (0 ≤ pred_rel_max ≤ 1)
        load["pred_rel_max"] = 0

        # Compensation for consuming less (i.e. voluntary demand reduction) (€/MWh)
        load["cost_red"] = 1000

        # Compensation for load curtailment (i.e. involuntary demand reduction) (€/MWh)
        load["cost_curt"] = 10000

        # Whether load is flexible (boolean)
        load["flex"] = 1

        # Power factor angle θ, giving the reactive power as Q = P ⨉ tan(θ)
        load["pf_angle"] = atan(load["qd"]/load["pd"])

        # Rescale cost and power input values to the p.u. values used internally in the model
        rescale_cost = x -> x*data["baseMVA"]
        _PM._apply_func!(load, "cost_red", rescale_cost)
        _PM._apply_func!(load, "cost_curt", rescale_cost)
    end
    return data
end

# This function selects the hours for the OPF calculation based on the input

function select_hours(year; selection = Dict{String, Any}("all" => true))
    if mod(year, 4) == 0
        all_hours = 1:(366 * 48)
    else
        all_hours = 1:(365 * 48)
    end
    
    if haskey(selection, "all") && selection["all"] == true
        selected_hours = all_hours
    elseif haskey(selection, "hour_range")
        selected_hours = all_hours[selection["hour_range"]]
    end

    return selected_hours
end
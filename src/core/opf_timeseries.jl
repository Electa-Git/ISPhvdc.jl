function prepare_mn_opf_data!(opf_mn_data, grid_data, demand_series, pv_series, wind_series, hours)
    for hour in hours
        opf_mn_data["nw"]["$hour"]["per_unit"] = true

        for (l, load) in opf_mn_data["nw"]["$hour"]["load"]
            load_bus = load["load_bus"]
            area_code = opf_mn_data["nw"]["$hour"]["bus"]["$load_bus"]["area"]
            area = opf_mn_data["nw"]["$hour"]["areas"]["$area_code"]
            demand_ratio = demand_series[area][hour]

            if load["pd"] >= 0 
                load["pd"] = grid_data["load"][l]["pd"] * demand_ratio * 0.9
                load["qd"] = grid_data["load"][l]["qd"] * demand_ratio * 0.9
            end
        end

        for (g, gen) in opf_mn_data["nw"]["$hour"]["gen"]
            trace = 1
            gen_bus = gen["gen_bus"]
            area_code = opf_mn_data["nw"]["$hour"]["bus"]["$gen_bus"]["area"]
            area = opf_mn_data["nw"]["$hour"]["areas"]["$area_code"]
            if gen["type"] == "Wind"
                trace = wind_series[area][hour] * 1.1
            elseif gen["type"] == "Solar"
                if !isempty(pv_series[area])
                    trace = pv_series[area][hour] * 1.1
                end
            end
            gen["pmax"] = grid_data["gen"][g]["pmax"] * trace
        end

        if haskey(opf_mn_data["nw"]["$hour"], "aggregated_data")
            delete!(opf_mn_data["nw"]["$hour"], "aggregated_data")
        end
    end
    
    return opf_mn_data
end


function prepare_hourly_opf_data!(hourly_data, grid_data, demand_series, pv_series, wind_series, rez_pv, rez_wind, time_stamp)

    hour = time_stamp

    hourly_data["pst"] = Dict{String, Any}()

    for (l, load) in hourly_data["load"]
        load_bus = load["load_bus"]
        area_code = hourly_data["bus"]["$load_bus"]["area"]
        area = hourly_data["areas"]["$area_code"]
        demand_ratio = demand_series[area][hour]

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
        load["cost_curt"] = 1e4 

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


function multi_network_uc_data(data, total_demand, total_dn_demand, pv_series, wind_series, rez_pv, rez_wind, hours, generator_contingencies; dn_res_factor = 0.5, no_dc_cont = false)
    data["pst"] = Dict{String, Any}()
    number_of_hours = length(hours)
    tie_line_contingencies = length(data["tie_lines"])
    conv_keys, dc_branch_keys = get_dc_continegncy_keys(data; no_dc_cont = no_dc_cont) 
    converter_contingencies = length(conv_keys) 
    dc_branch_contingencies = length(dc_branch_keys) 
    number_of_contingencies = sum(generator_contingencies) + tie_line_contingencies + converter_contingencies +  dc_branch_contingencies + 1 # to also add the N case

    # This for loop determines which "network" belongs to an hour, and which to a contingency, for book-keeping of the network ids
    # Format: [h1, c1 ... cn, h2, c1 ... cn, .... , hn, c1 ... cn]
    hour_ids = [];
    cont_ids = [];
    for i in 1:number_of_hours * number_of_contingencies
        if mod(i, number_of_contingencies) == 1
            push!(hour_ids, i)
        else
            push!(cont_ids, i)
        end
    end
    mn_data = _IM.replicate(data, number_of_hours * number_of_contingencies, Set{String}(["source_type", "name", "source_version", "per_unit"]))
    mn_data["hour_ids"] = hour_ids
    mn_data["cont_ids"] = cont_ids
    mn_data["number_of_hours"] = number_of_hours
    mn_data["number_of_contingencies"] = number_of_contingencies

    for idx in 1:number_of_hours
        hour = hours[idx]
        nw_start = 1 + (idx - 1) * (number_of_contingencies)
        nw_ids = nw_start:(nw_start + number_of_contingencies-1)
        for nw in nw_ids
            for (l, load) in mn_data["nw"]["$nw"]["load"]
                load_bus = load["load_bus"]
                area_code = data["bus"]["$load_bus"]["area"]
                area = data["areas"]["$area_code"]
                wind_dn = total_dn_demand["Wind"][area][hour]
                if haskey(total_dn_demand["PV"], area)
                    pv_dn = total_dn_demand["PV"][area][hour]
                else
                    pv_dn = 0
                end
                demand_trace = total_demand[area][hour] + dn_res_factor * (pv_dn + wind_dn)

                if load["pd"] >= 0 
                    q_ratio = data["load"][l]["qd"] / data["load"][l]["pd"]
                    load["pd"] = data["load"][l]["pd"] / data["aggregated_data"][area]["demand"] * demand_trace
                    load["qd"] = data["load"][l]["pd"] * q_ratio
                end
            end
        end
    end
    for idx in 1:number_of_hours
        hour = hours[idx]
        nw_start = 1 + (idx - 1) * (number_of_contingencies)
        nw_ids = nw_start:(nw_start + number_of_contingencies-1)
        for nw in nw_ids
            for (g, gen) in mn_data["nw"]["$nw"]["gen"]
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
                    area_code = data["bus"]["$gen_bus"]["area"]
                    area = data["areas"]["$area_code"]
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
                gen["pmax"] = data["gen"][g]["pmax"] * trace
            end
        end
    end

    create_contingencies!(mn_data, generator_contingencies, conv_keys, dc_branch_keys)

    for (n, network) in mn_data["nw"]
        network["zones"] = Dict{String, Any}("1" => Dict("source_id" => Any["zones", 1], "zone" => 1, "index" => 1), "2" => Dict("source_id" => Any["zones", 2], "zone" => 5, "index" => 2))
        network["areas"] = Dict{String, Any}("1" => Dict("source_id" => Any["areas", 1], "area" => 1, "index" => 1), "2" => Dict("source_id" => Any["areas", 2], "area" => 3, "index" => 2), "3" => Dict("source_id" => Any["areas", 3], "area" => 4, "index" => 3), "4" => Dict("source_id" => Any["areas", 4], "area" => 5, "index" => 4))
        for (b, bus) in network["bus"]
            if bus["area"] == 5
                bus["zone"] = 5
            elseif bus["area"] == 2
                bus["area"] = 1
            end
        end
        for (c, conv) in network["convdc"]
            conv_bus = conv["busac_i"]
            conv["zone"] = network["bus"]["$conv_bus"]["zone"]
            conv["area"] = network["bus"]["$conv_bus"]["area"]
        end
        for (g, gen) in network["gen"]
            if gen["gen_status"] == 1
                gen_bus = gen["gen_bus"]
                gen["zone"] = network["bus"]["$gen_bus"]["zone"]
                gen["area"] = network["bus"]["$gen_bus"]["area"]
            end
        end
    end

    for (n, network) in mn_data["nw"]
        if haskey(network, "aggregated_data")
            delete!(network, "aggregated_data")
        end
    end
    
    return mn_data
end

function create_contingencies!(mn_data, generator_contingencies, conv_keys, dc_branch_keys)
    for hour in mn_data["hour_ids"]
        gen_keys = find_most_severe_generator_contingecnies(mn_data, hour, generator_contingencies)
        tie_line_keys = sort(parse.(Int, collect(keys(mn_data["nw"]["1"]["tie_lines"]))))

        mn_data["nw"]["$hour"]["contingency"] = Dict{String, Any}("gen_id" => nothing, "branch_id" => nothing, "conv_id" => nothing, "dcbranch_id" => nothing)
        for idx in 1:length(gen_keys)
            nw = hour + idx
            gen_id = gen_keys[idx]
            mn_data["nw"]["$nw"]["contingency"] = Dict{String, Any}("gen_id" => gen_id, "branch_id" => nothing, "conv_id" => nothing, "dcbranch_id" => nothing)
        end
        for idx in 1:length(tie_line_keys)
            nw = hour + length(gen_keys) +idx
            line_id = tie_line_keys[idx]
            mn_data["nw"]["$nw"]["contingency"] = Dict{String, Any}("gen_id" => nothing, "branch_id" => line_id, "conv_id" => nothing, "dcbranch_id" => nothing)
        end
        for idx in 1:length(conv_keys)
            nw = hour + length(gen_keys) + length(tie_line_keys) +idx
            conv_id = conv_keys[idx]
            mn_data["nw"]["$nw"]["contingency"] = Dict{String, Any}("gen_id" => nothing, "branch_id" => nothing, "conv_id" => conv_id, "dcbranch_id" => nothing)
        end
        for idx in 1:length(dc_branch_keys)
            nw =hour + length(gen_keys) + length(tie_line_keys) + length(conv_keys) +idx
            line_id = dc_branch_keys[idx]
            mn_data["nw"]["$nw"]["contingency"] = Dict{String, Any}("gen_id" => nothing, "branch_id" => nothing, "conv_id" => nothing, "dcbranch_id" => line_id)
        end
    end
    return mn_data
end

function find_most_severe_generator_contingecnies(mn_data, hour, generator_contingencies)

    area_gens = Dict{String, Any}(["$area" => [] for area in 1:5])
    area_gen_keys = []

    for (g, gen) in mn_data["nw"]["$hour"]["gen"]
        area = gen["area"]
        push!(area_gens["$area"], g)
    end

    for (area, gens) in area_gens
        gen_keys = sort(parse.(Int, collect(gens)))
        pmax = [mn_data["nw"]["$hour"]["gen"]["$g"]["pmax"] for g in gen_keys]
        if !isempty(pmax)
            contingencies =  generator_contingencies[parse(Int, area)]
            keys = gen_keys[sortperm(pmax, rev = true)[1:contingencies]]
            for key in keys
              push!(area_gen_keys, key)
            end
        end
    end

    return area_gen_keys
end


function multi_network_capacity_data(data, total_demand, total_dn_demand, pv_series, wind_series, rez_pv, rez_wind, hours; dn_res_factor = 0.5)
    data["pst"] = Dict{String, Any}()
    number_of_hours = length(hours)

    mn_data = _IM.replicate(data, number_of_hours, Set{String}(["source_type", "name", "source_version", "per_unit"]))

    for nw in 1:length(hours)
        hour = hours[nw]
        for (l, load) in mn_data["nw"]["$nw"]["load"]
            load_bus = load["load_bus"]
            area_code = data["bus"]["$load_bus"]["area"]
            area = data["areas"]["$area_code"]
            wind_dn = total_dn_demand["Wind"][area][hour]
            if haskey(total_dn_demand["PV"], area)
                pv_dn = total_dn_demand["PV"][area][hour]
            else
                pv_dn = 0
            end
            demand_trace = total_demand[area][hour] + dn_res_factor * (pv_dn + wind_dn)

            if load["pd"] >= 0 
                q_ratio = data["load"][l]["qd"] / data["load"][l]["pd"]
                load["pd"] = data["load"][l]["pd"] / data["aggregated_data"][area]["demand"] * demand_trace
                load["qd"] = data["load"][l]["pd"] * q_ratio
            end
        end
    end
    for nw in 1:length(hours)
        hour = hours[nw]
        for (g, gen) in mn_data["nw"]["$nw"]["gen"]
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
                area_code = data["bus"]["$gen_bus"]["area"]
                area = data["areas"]["$area_code"]
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
            gen["pmax"] = data["gen"][g]["pmax"] * trace
        end
    end

    for (n, network) in mn_data["nw"]
        if haskey(network, "aggregated_data")
            delete!(network, "aggregated_data")
        end
    end
        
    return mn_data
end

function get_dc_continegncy_keys(data; no_dc_cont = false)
    if no_dc_cont
        conv_keys = []
        dc_branch_keys = []
    else
        conv_keys = []
        for (c, conv) in data["convdc"]
            area_1, area_2 = find_areas_dc(data, conv)
            if haskey(conv, "rez_connection") && conv["rez_connection"] == true
                push!(conv_keys, parse(Int, c))
            elseif area_1 !== area_2
                push!(conv_keys, parse(Int, c))
            end
        end
    end
        conv_keys = sort(conv_keys)
        dc_branch_keys = []
    #     for (b, br) in data["branchdc"]
    #         # if !haskey(br, "rez_connection")
    #             push!(dc_branch_keys, parse(Int, b))
    #         # elseif  haskey(br, "rez_connection") && br["rez_connection"] == false
    #         #     push!(dc_branch_keys, parse(Int, b))
    #         # end
    #     end
    #     dc_branch_keys = sort(dc_branch_keys)
    # end
    return conv_keys, dc_branch_keys
end

function find_areas_dc(data, conv)
    ac_bus = conv["busac_i"]
    area_1 = data["bus"]["$ac_bus"]["area"]
    dc_bus = conv["busdc_i"]
    dc_bus_ = []
    for (b, br) in data["branchdc"]
        if br["fbusdc"] == dc_bus
            dc_bus_ = br["tbusdc"]
        elseif br["tbusdc"] == dc_bus
            dc_bus_ = br["fbusdc"]
        end 
    end
    ac_bus_ = []
    for (c, conv) in data["convdc"]
        if conv["busdc_i"] == dc_bus_
            ac_bus_ = conv["busac_i"]
        end
    end
    area_2 = data["bus"]["$ac_bus_"]["area"]

    if area_1 == 2
        area_1 = 1
    end
    if area_2 == 2
        area_2 = 1
    end
    
    return area_1, area_2
end
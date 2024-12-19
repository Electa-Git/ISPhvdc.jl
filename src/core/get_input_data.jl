function download_isp_data()
    # Download ISP input workbook
    mkdir(joinpath("data", "inputs"))
    filename_inputs = joinpath("data", "Inputs", "Inputs assumptions and scenarios workbook.xlsx")
    url_inputs = "https://aemo.com.au/-/media/files/major-publications/isp/2022/2022-documents/inputs-assumptions-and-scenarios-workbook.xlsx?la=en"
    Downloads.download(url_inputs, filename_inputs)


    filename_solar = joinpath("data", "2022-isp-solar-traces.zip")
    url_solar = "https://aemo.com.au/-/media/files/major-publications/isp/2022/2022-documents/2022-isp-solar-traces.zip?la=en"
    Downloads.download(url_solar, filename_solar)
    InfoZIP.unzip(filename_solar, "data")
    rm(filename_solar)

    filename_wind = joinpath("data", "2022 ISP Wind traces.zip")
    url_wind = "https://aemo.com.au/-/media/files/major-publications/isp/2022/2022-documents/2022-isp-wind-traces.zip?la=en"
    Downloads.download(url_wind, filename_wind)
    InfoZIP.unzip(filename_wind, "data")
    rm(filename_wind)

    filename_generation = joinpath("data", "Generation outlook.zip")
    url_generation = "https://aemo.com.au/-/media/files/major-publications/isp/2022/2022-documents/generation-outlook.zip?la=en"
    Downloads.download(url_generation, filename_generation)
    mkdir(joinpath("data", "Generation outlook"))
    InfoZIP.unzip(filename_generation, "data/Generation outlook")
    rm(filename_generation)

    filename_model = joinpath("data", "2022 ISP model.zip")
    url_model = "https://aemo.com.au/-/media/files/major-publications/isp/2022/2022-documents/2022-isp-model.zip?la=en"
    Downloads.download(url_model, filename_model)
    InfoZIP.unzip(filename_model, "data")
    rm(filename_model)
end

function get_demand_data(scenario, year; data_dir="data")
    data_folder = joinpath(data_dir, "2022 Final ISP Model", scenario, "Traces", "demand")
    all_demand_files = readdir(data_folder)

    demand = Dict{String,Any}()
    demand_series = Dict{String,Any}()

    for file in all_demand_files
        region = file[1:(collect(findfirst("_", file)).-1)[1]]
        print("Reading demand trace for ", region, "\n")
        demand_region = _DF.DataFrame(CSV.File(joinpath(data_folder, file)))
        demand[region] = demand_region[demand_region[!, :Year].==year, :]
    end

    demand_states = aggregate_demand(demand)


    for (d, demand) in demand_states
        demand_series[d] = zeros(1, 48 * _DF.nrow(demand))
        for row in 1:_DF.nrow(demand)
            demand_series[d][(((row-1)*48)+1):(row*48)] = collect(values(demand[row, 4:end]))
        end
    end

    return demand_series
end


function get_dn_demand_data(scenario, year; data_dir="data")
    data_folder = joinpath(data_dir, "2022 Final ISP Model", scenario, "Traces", "load_subtractor")
    all_demand_files = readdir(data_folder)

    demand = Dict{String,Any}("PV" => Dict{String,Any}(), "Wind" => Dict{String,Any}())

    for file in all_demand_files
        region = file[(collect(findall("_", file)[2])[1]+1):(collect(findall("_", file)[3])[1]-1)]
        type = file[(collect(findall("_", file)[1])[1]+1):(collect(findall("_", file)[2])[1]-1)]
        print("Reading ", type, " load subtractor trace for ", region, "\n")
        demand_region = _DF.DataFrame(CSV.File(joinpath(data_folder, file)))
        demand[type][region] = demand_region[demand_region[!, :Year].==year, :]
    end

    demand_series = Dict{String,Any}()
    for (t, type) in demand
        demand_series[t] = Dict{String,Any}()
        for (s, state) in type
            demand_series[t][s] = zeros(1, 48 * _DF.nrow(state))
            for row in 1:_DF.nrow(state)
                demand_series[t][s][(((row-1)*48)+1):(row*48)] = collect(values(state[row, 4:end]))
            end
        end
    end

    return demand_series
end

function aggregate_demand(demand)

    aggregated_demand = Dict{String,Any}("NSW" => demand["NNSW"], "VIC" => demand["VIC"], "QLD" => demand["SQ"], "SA" => demand["SA"], "TAS" => demand["TAS"])

    for i in 1:size(aggregated_demand["NSW"], 1)
        for j in 4:size(aggregated_demand["NSW"], 2)
            aggregated_demand["NSW"][i, j] = demand["NNSW"][i, j] + demand["SNSW"][i, j] + demand["CNSW"][i, j] + demand["SNW"][i, j]
            aggregated_demand["QLD"][i, j] = demand["SQ"][i, j] + demand["CNQ"][i, j] + demand["GG"][i, j]
        end
    end

    return aggregated_demand
end

function get_rez_capacity_data(scenario, year, data_dir; CDP="CDP10")
    data_folder = joinpath(data_dir, "Generation Outlook", "Final ISP Results", "Scenarios")
    file_name = joinpath(data_folder, join(["2022 Final ISP results workbook - ", scenario[10:end], " - Updated Inputs.xlsx"]))
    rez_all_data = XLSX.readtable(file_name, "REZ Generation Capacity", "A:AG", first_row=3, header=true, stop_in_empty_row=false,) |> _DF.DataFrame # read all data
    rez_all_data = nothing
    try
        rez_all_data = _DF.DataFrame(
            XLSX.readtable(file_name, "REZ Generation Capacity", "A:AG", first_row=3, header=true, stop_in_empty_row=false,)...
        ) # when called with v0.7.10
    catch
        rez_all_data = _DF.DataFrame(
            XLSX.readtable(file_name, "REZ Generation Capacity", "A:AG", first_row=3, header=true, stop_in_empty_row=false,)
        ) # when called with v10.2
    end
    filter!(row -> isequal(row.CDP, CDP), rez_all_data) # filter for the CDP

    # extract the data for solar, onshore wind and offshore wind
    rez_pv = Matrix(filter(row -> row.Technology == "Solar", rez_all_data))
    rez_onshore_wind = Matrix(filter(row -> row.Technology == "Wind", rez_all_data))
    rez_offshore_wind = Matrix(filter(row -> row.Technology == "Offshore Wind", rez_all_data))

    REZ_capacity = Dict{String,Any}()

    year_idx = (year - 2024) + 6

    if year > 2051
        print("The selected year is out of range, installed capacities are provided between 2024 and 2051", "\n")
    else
        REZ_capacity["pv"] = rez_pv[:, [1, 2, 3, 4, 5, year_idx]]
        REZ_capacity["onshore_wind"] = rez_onshore_wind[:, [1, 2, 3, 4, 5, year_idx]]
        REZ_capacity["offshore_wind"] = rez_offshore_wind[:, [1, 2, 3, 4, 5, year_idx]]
    end

    return REZ_capacity
end


function get_rez_grid_extensions(; data_dir="data")
    file_name = joinpath(data_dir, "NEM_REZ_extensions.xlsx")
    ac = XLSX.readdata(file_name, "AC", "A1:S34")
    dc = XLSX.readdata(file_name, "DC", "A1:P5")

    rez_extensions = Dict{String,Any}()
    rez_extensions["ac"] = ac
    rez_extensions["dc"] = dc
    return rez_extensions
end

function get_generator_information(; data_dir="data")
    file_name = joinpath(data_dir, "inputs", "Inputs assumptions and scenarios workbook.xlsx")
    gen_info_existing = XLSX.readdata(file_name, "Existing Gen Data Summary", "B12:U229")
    gen_info_committed = XLSX.readdata(file_name, "Existing Gen Data Summary", "B232:U277")
    gen_info_anticipated = XLSX.readdata(file_name, "Existing Gen Data Summary", "B281:U294")

    gen_info = [gen_info_existing; gen_info_committed; gen_info_anticipated]

    remove_extra_info!(gen_info)

    return gen_info
end

function remove_extra_info!(gen_info)
    for idx = 1:size(gen_info, 1)
        gen_name = gen_info[idx, 1]
        if !isnothing(findfirst("Solar", gen_name))
            solar_idx = collect(findfirst("Solar", gen_name))[1] - 2
            gen_info[idx, 1] = gen_name[1:solar_idx]
        elseif !isnothing(findfirst("Wind", gen_name))
            wind_idx = collect(findfirst("Wind", gen_name))[1] - 2
            gen_info[idx, 1] = gen_name[1:wind_idx]
        end
    end

    return gen_info
end

function get_res_timeseries(year; data_dir="data")
    pv = get_pv_timeseries(year; data_dir)
    wind = get_wind_timeseries(year; data_dir)

    return pv, wind
end


function get_pv_timeseries(year; data_dir="data")
    data_folder = joinpath(data_dir, "solar")
    all_files = readdir(data_folder)

    pv = Dict{String,Any}()

    for file in all_files

        if !isempty(findall("SAT", file))
            plant = file[1:(findall("SAT", file)[1][1]-2)]
        elseif !isempty(findall("CST", file))
            plant = file[1:(findall("CST", file)[1][1]-2)]
        elseif !isempty(findall("FFP", file))
            plant = file[1:(findall("FFP", file)[1][1]-2)]
        elseif !isempty(findall("PEG", file))
            plant = file[1:(findall("PEG", file)[1][1]-2)]
        end
        print(file)
        print("Reading PV trace for ", plant, "\n")
        pv_plant = _DF.DataFrame(CSV.File(joinpath(data_folder, file)))
        pv[plant] = pv_plant[pv_plant[!, :Year].==year, :]
    end

    return pv
end

function get_wind_timeseries(year; data_dir="data")
    data_folder = joinpath(data_dir, "wind")
    all_files = readdir(data_folder)

    wind = Dict{String,Any}()

    for file in all_files
        plant = file[1:(findall("Ref", file)[1][1]-2)]
        print("Reading Wind trace for ", plant, "\n")
        wind_plant = _DF.DataFrame(CSV.File(joinpath(data_folder, file)))
        wind[plant] = wind_plant[wind_plant[!, :Year].==year, :]
    end

    return wind
end


function aggregate_res_timeseries(res_ts, generator_info, type::String)

    aggregated_res = Dict{String,Any}("NSW" => Dict{String,Any}("values" => _DF.DataFrame(), "count" => 0), "VIC" => Dict{String,Any}("values" => _DF.DataFrame(), "count" => 0),
        "QLD" => Dict{String,Any}("values" => _DF.DataFrame(), "count" => 0), "SA" => Dict{String,Any}("values" => _DF.DataFrame(), "count" => 0), "TAS" => Dict{String,Any}("values" => _DF.DataFrame(), "count" => 0))

    res_series = Dict{String,Any}("NSW" => [], "VIC" => [], "QLD" => [], "SA" => [], "TAS" => [])
    count = Dict{String,Any}("NSW" => 0, "VIC" => 0, "QLD" => 0, "SA" => 0, "TAS" => 0)
    for (plant, profile) in res_ts
        if plant[1:3] !== "REZ" && plant[3] !== '_'
            similarity = [_SD.compare(lowercase(plant), lowercase(generator_info[idx, 1]), _SD.Jaro()) for idx in findall(generator_info[:, 6] .== type)]
            ids = [idx for idx in findall(generator_info[:, 6] .== type)]
            best_fit = findmax(similarity)[2][1]
            name = generator_info[ids[best_fit], 1]
            state = generator_info[ids[best_fit], 3]
            print(plant, " -> ", name, "\n")
            count[state] = count[state] + 1
            if isempty(res_series[state])
                res_series[state] = zeros(1, 48 * _DF.nrow(profile))
                for row in 1:_DF.nrow(profile)
                    res_series[state][(((row-1)*48)+1):(row*48)] = collect(values(profile[row, 4:end]))
                end
            else
                for row in 1:_DF.nrow(profile)
                    res_series[state][(((row-1)*48)+1):(row*48)] = res_series[state][(((row-1)*48)+1):(row*48)] .+ collect(values(profile[row, 4:end]))
                end
            end
        end
    end

    for (s, state) in res_series
        res_series[s] = res_series[s] ./ count[s]
    end

    return res_series, count

end

function make_rez_time_series(res_ts)
    ts = Dict{String,Any}()
    rez_name = nothing
    for (plant, profile) in res_ts
        if plant[1:3] == "REZ"
            rez_name = plant[5:6]
        elseif plant[3] == '_'
            rez_name = plant[1:2]
        end
        if !isnothing(rez_name)
            res_series = zeros(1, 48 * _DF.nrow(profile))
            for row in 1:_DF.nrow(profile)
                res_series[(((row-1)*48)+1):(row*48)] = collect(values(profile[row, 4:end]))
            end
            push!(ts, rez_name => res_series)
        end
    end
    return ts
end

function add_rez_and_connections!(data, extensions, rez)

    for row in eachrow(extensions["ac"])
        if row[1] !== "REZ ID"
            if row[15] !== 0
                tbus = row[15]
            else
                tbus = maximum(sort(parse.(Int, keys(data["bus"])))) + 1
                add_ac_bus!(data, row, tbus)
            end
            add_ac_branch!(data, row, tbus)
            add_generator!(data, row, tbus, rez)
        end
    end

    for row in eachrow(extensions["dc"])
        if row[1] !== "REZ ID"
            fbus = row[11]
            if row[12] !== 0
                tbus = row[12]
            else
                tbus = maximum(sort(parse.(Int, keys(data["bus"])))) + 1
                add_ac_bus!(data, row, tbus)
            end
            fbusdc = maximum(sort(parse.(Int, keys(data["busdc"])))) + 1
            tbusdc = maximum(sort(parse.(Int, keys(data["busdc"])))) + 2
            add_dc_bus!(data, row, fbusdc)
            add_dc_bus!(data, row, tbusdc)
            if row[1] !== "MARINUS"
                add_dc_converter!(data, row, fbus, fbusdc; rez_connection=true)
                add_dc_converter!(data, row, tbus, tbusdc; rez_connection=true)
                add_dc_branch!(data, row, fbusdc, tbusdc; rez_connection=true)
            else
                add_dc_converter!(data, row, fbus, fbusdc; rez_connection=false)
                add_dc_converter!(data, row, tbus, tbusdc; rez_connection=false)
                add_dc_branch!(data, row, fbusdc, tbusdc; rez_connection=false)
            end
            if row[1] !== "MARINUS"
                add_generator!(data, row, tbus, rez)
            end
        end
    end

    return data
end


function add_ac_branch!(data, row, tbus)
    for circuits in 1:row[8]
        branch_id = maximum(sort(parse.(Int, keys(data["branch"])))) + 1
        branch = Dict{String,Any}()
        fbus = row[14]
        branch["f_bus"] = fbus
        branch["t_bus"] = tbus
        base_kv = row[6]
        base_mva = data["baseMVA"]
        branch["br_r"] = to_pu_z(row[11], base_kv, base_mva)
        branch["br_x"] = to_pu_z(row[12], base_kv, base_mva)
        branch["b_fr"] = to_pu_y(row[13], base_kv, base_mva) / 2
        branch["b_to"] = to_pu_y(row[13], base_kv, base_mva) / 2
        branch["g_fr"] = branch["g_to"] = 0.0
        branch["br_status"] = 1
        branch["rate_a"] = branch["rate_b"] = branch["rate_c"] = row[7] / base_mva
        branch["name"] = join(["connection ", row[1, 1]])
        branch["index"] = branch_id
        branch["angmin"] = -pi / 2
        branch["angmax"] = pi / 2
        branch["transformer"] = row[10]
        branch["tap"] = 1.0
        branch["shift"] = 0.0
        push!(data["branch"], "$branch_id" => branch)
    end
end

function to_pu_z(z, base_kv, base_mva)
    z_base = (base_kv * 1e3)^2 / (base_mva * 1e6)
    z_pu = z / z_base
    return z_pu
end

function to_pu_y(y, base_kv, base_mva)
    y_base = (base_mva * 1e6) / (base_kv * 1e3)^2
    y_pu = y / y_base
    return y_pu
end

function add_ac_bus!(data, row, tbus)
    bus = Dict{String,Any}()
    bus["zone"] = 1
    bus["index"] = tbus
    bus["bus_i"] = tbus
    bus["bus_type"] = 2
    bus["name"] = join(["connection ", row[1, 1]])
    bus["va"] = 0.0
    bus["vm"] = 1.0
    bus["vmax"] = 1.1
    bus["vmin"] = 0.9
    bus["base_kv"] = row[6]
    if length(row) > 16
        bus["lat"] = row[18]
        bus["lon"] = row[19]
    else
        bus["lat"] = row[15]
        bus["lat"] = row[16]
    end
    area_ = 0
    for (area_idx, area) in data["areas"]
        if area[1] == row[1][1]
            area_ = parse(Int, area_idx)
        end
    end
    bus["area"] = area_

    push!(data["bus"], "$tbus" => bus)
end

function add_generator!(data, row, tbus, rez_capacities; max_gen_power=1000)
    if any(row[1] .== rez_capacities["pv"][:, 3])
        rez_id = findfirst(row[1] .== rez_capacities["pv"][:, 3])
        rez_power = rez_capacities["pv"][rez_id, :][end]
        gens = Int(floor(rez_power / max_gen_power))
        for idx in 1:(gens+1)
            gen_id = maximum(sort(parse.(Int, keys(data["gen"])))) + 1
            gen_power = max(0, min(max_gen_power, (rez_power - (idx * max_gen_power))))
            gen = generator_data(tbus, rez_capacities["pv"][rez_id, :], gen_power, gen_id, "Solar", data["baseMVA"])
            push!(data["gen"], "$gen_id" => gen)
        end
    end
    if any(row[1] .== rez_capacities["onshore_wind"][:, 3])
        rez_id = findfirst(row[1] .== rez_capacities["onshore_wind"][:, 3])
        rez_power = rez_capacities["onshore_wind"][rez_id, :][end]
        gens = Int(floor(rez_power / max_gen_power))
        for idx in 1:(gens+1)
            gen_id = maximum(sort(parse.(Int, keys(data["gen"])))) + 1
            gen_power = max(0, min(max_gen_power, (rez_power - (idx * max_gen_power))))
            gen = generator_data(tbus, rez_capacities["onshore_wind"][rez_id, :], gen_power, gen_id, "Wind", data["baseMVA"])
            push!(data["gen"], "$gen_id" => gen)
        end
    end

    # To Do: Add offshore wind
end

function generator_data(tbus, rez, power, gen_id, type, basemva)
    gen = Dict{String,Any}()
    gen["mbase"] = gen["pmax"] = power / basemva
    gen["qmax"] = power / 2 / basemva
    gen["qmin"] = -power / 2 / basemva
    gen["name"] = rez[3]
    gen["pmin"] = gen["qg"] = gen["pg"] = 0
    gen["ncost"] = 2
    gen["model"] = 2
    gen["gen_bus"] = tbus
    gen["index"] = gen_id
    gen["fuel"] = type
    gen["cost"] = [1.0 0.0] # check for correct costs later
    gen["gen_status"] = 1
    gen["type"] = type

    return gen
end


function add_dc_converter!(data, row, acbus, dcbus_id; rez_connection=false)
    for circuits in 1:row[8]
        conv_id = maximum(sort(parse.(Int, keys(data["convdc"])))) + 1
        conv = Dict{String,Any}()

        conv["busac_i"] = acbus
        conv["busdc_i"] = dcbus_id
        conv["type_dc"] = 2
        conv["type_ac"] = 1
        conv["index"] = conv_id
        conv["P_g"] = conv["Q_g"] = conv["dVdcset"] = conv["Pdcset"] = conv["droop"] = 0
        conv["islcc"] = 0
        conv["status"] = 1
        conv["Vdcset"] = conv["Vtar"] = conv["tm"] = 1
        conv["filter"] = conv["reactor"] = conv["transformer"] = 0
        conv["Vmmax"] = 1.1
        conv["Vmmin"] = 0.9

        conv["Pacmax"] = row[7] / data["baseMVA"]
        conv["Pacmin"] = -row[7] / data["baseMVA"]
        conv["Qacmax"] = row[7] / 2 / data["baseMVA"]
        conv["Qacmin"] = -row[7] / 2 / data["baseMVA"]
        conv["Pacmax"] = row[7] / data["baseMVA"]
        conv["Pacrated"] = conv["Imax"] = 1.1 * row[7] / data["baseMVA"]
        conv["Qacrated"] = 1.1 * row[7] / 2 / data["baseMVA"]
        conv["basekVac"] = row[6]

        acbus_voltage = data["bus"]["$acbus"]["base_kv"]
        zbase = (acbus_voltage * 1e3)^2 / (data["baseMVA"] * 1e6)
        zbase_dc = (row[6] * 1e3)^2 / (data["baseMVA"] * 1e6)
        xr_trafo = 35
        sc_trafo = 0.15
        xr_reactor = 30
        sc_reactor = 7.5

        ztf = (acbus_voltage * 1e3)^2 / (row[7] * 1e6) * sc_trafo
        zreactor = (acbus_voltage * 1e3)^2 / (row[7] * 1e6) * sc_reactor
        conv["rtf"] = ztf * cos(atan(xr_trafo)) / zbase
        conv["xtf"] = ztf * sin(atan(xr_trafo)) / zbase
        conv["rc"] = zreactor * cos(atan(xr_reactor)) / zbase
        conv["xc"] = zreactor * sin(atan(xr_reactor)) / zbase
        conv["bf"] = 0
        conv["LossA"] = 1.1033 * 1e-3 * row[7] / data["baseMVA"]
        conv["LossB"] = 0.0035
        conv["LossCinv"] = conv["LossCrec"] = 0.0035 / row[7] * zbase_dc

        conv["rez_connection"] = rez_connection # used for removing these from the set of contingencies.....


        push!(data["convdc"], "$conv_id" => conv)
    end

end

function add_dc_branch!(data, row, dcbus_id1, dcbus_id2; rez_connection=false)
    for circuits in 1:row[8]
        dcbranch_id = maximum(sort(parse.(Int, keys(data["branchdc"])))) + 1
        dcbranch = Dict{String,Any}()
        dcbranch["fbusdc"] = dcbus_id1
        dcbranch["tbusdc"] = dcbus_id2
        dcbranch["rateA"] = dcbranch["rateB"] = dcbranch["rateC"] = row[7] / data["baseMVA"]
        dcbranch["index"] = dcbranch_id
        dcbranch["l"] = dcbranch["c"] = 0
        dcbranch["status"] = 1
        dcvoltage = data["busdc"]["$dcbus_id1"]["basekVdc"]
        zbase = (dcvoltage * 1e3)^2 / (data["baseMVA"] * 1e6)
        dcbranch["r"] = row[10] / zbase
        dcbranch["rez_connection"] = rez_connection

        push!(data["branchdc"], "$dcbranch_id" => dcbranch)
    end
end

function add_dc_bus!(data, row, dcbus)
    busdc = Dict{String,Any}()
    busdc["basekVdc"] = row[6]
    busdc["index"] = busdc["busdc_i"] = dcbus
    busdc["Vdc"] = 1
    busdc["Cdc"] = busdc["Pdc"] = 0
    busdc["grid"] = 1
    busdc["Vdcmax"] = 1.1
    busdc["Vdcmin"] = 0.9

    push!(data["busdc"], "$dcbus" => busdc)
end

function fix_data!(data)

    # Find isolated buses and put their demand zero
    bus_arcs = Dict{String,Any}([b => [] for (b, bus) in data["bus"]])
    for (b, branch) in data["branch"]
        fbus = branch["f_bus"]
        tbus = branch["t_bus"]
        if haskey(bus_arcs, "$fbus")
            push!(bus_arcs["$fbus"], parse(Int, b))
        end
        if haskey(bus_arcs, "$tbus")
            push!(bus_arcs["$tbus"], parse(Int, b))
        end
        branch["tap"] = 1.0
        branch["shift"] = 0.0
    end
    for (c, conv) in data["convdc"]
        cbus = conv["busac_i"]
        push!(bus_arcs["$cbus"], parse(Int, c))
    end

    for (l, load) in data["load"]
        load_bus = load["load_bus"]
        if isempty(bus_arcs["$load_bus"])
            load["pd"] = 0.0
            load["qd"] = 0.0
        end
    end

    # generator data comming from matlab model seems two orders of magnitude too small
    for (g, gen) in data["gen"]
        gen["cost"] = gen["cost"] .* data["baseMVA"]
    end

    return data
end

# function fix_hvdc_data_issues!(data; no_bass = false, no_terra = false, no_murray = false)

#     if no_bass == false
#         # BASS LINK 
#         data["branch"]["543"]["br_status"] = 0
#         delete!(data["bus"], "2113")
#         data["gen"]["264"]["gen_status"] = 0
#         data["gen"]["265"]["gen_status"] = 0
#     end

#     if no_terra == false
#         # TERRANORALINK
#         data["branch"]["1568"]["br_status"] = 0
#         data["branch"]["1569"]["br_status"] = 0
#         data["branch"]["1570"]["br_status"] = 0
#     end

#     if no_murray ==false
#         # MURRAYLINK
#         data["branch"]["1949"]["br_status"] = 0
#         data["gen"]["222"]["gen_status"] = 0
#     end

#     return data
# end

function fix_hvdc_data_issues!(data; no_bass=false, no_murray=false, no_terra=false)

    if no_bass == false
        # BASS LINK
        delete!(data["bus"], "2113")
        delete!(data["branch"], "543")
        delete!(data["shunt"], "237")
        delete!(data["gen"], "264")
        delete!(data["gen"], "265")
        data["bus"]["2250"]["bus_type"] = 1
    end

    if no_murray == false
        # MURRAYLINK
        delete!(data["gen"], "222")
        data["bus"]["986"]["bus_type"] = 1
        delete!(data["load"], "395")
        delete!(data["shunt"], "98")
    end

    if no_terra == false
        # TERRANORALINK
        delete!(data["gen"], "206")
        data["bus"]["211"]["bus_type"] = 1
        delete!(data["load"], "75")
    end
    return data
end











function add_area_dict!(data_nem)
    data_nem["areas"] = Dict{String,Any}()
    data_nem["areas"]["1"] = "NSW"
    data_nem["areas"]["2"] = "VIC"
    data_nem["areas"]["3"] = "QLD"
    data_nem["areas"]["4"] = "SA"
    data_nem["areas"]["5"] = "TAS"
end


function merge_parallel_lines(data)
    ft_buses = Dict{String,Any}()
    for (b, branch) in data["branch"]
        ftbus = join([branch["f_bus"], "-", branch["t_bus"]])
        tfbus = join([branch["t_bus"], "-", branch["f_bus"]])
        if !haskey(ft_buses, ftbus) && !haskey(ft_buses, tfbus)
            push!(ft_buses, ftbus => b)
        elseif haskey(ft_buses, ftbus)
            branch_id = ft_buses[ftbus]
            merge_branches!(data, branch_id, branch)
            delete!(data["branch"], b)
        elseif haskey(ft_buses, tfbus)
            branch_id = ft_buses[tfbus]
            merge_branches!(data, branch_id, branch)
            delete!(data["branch"], b)
        end
    end

    return data
end

function merge_branches!(data, b_idx, parallel_br)
    data["branch"][b_idx]["rate_a"] = data["branch"][b_idx]["rate_a"] + parallel_br["rate_a"]
    data["branch"][b_idx]["rate_b"] = data["branch"][b_idx]["rate_b"] + parallel_br["rate_b"]
    data["branch"][b_idx]["rate_c"] = data["branch"][b_idx]["rate_c"] + parallel_br["rate_c"]
    if haskey(data["branch"][b_idx], "c_rating") && haskey(parallel_br, "c_rating")
        data["branch"][b_idx]["c_rating"] = data["branch"][b_idx]["c_rating"] + parallel_br["c_rating"]
    end
    data["branch"][b_idx]["br_r"] = (data["branch"][b_idx]["br_r"] * parallel_br["br_r"]) / (data["branch"][b_idx]["br_r"] + parallel_br["br_r"])
    data["branch"][b_idx]["br_x"] = (data["branch"][b_idx]["br_x"] * parallel_br["br_x"]) / (data["branch"][b_idx]["br_x"] + parallel_br["br_x"])
    data["branch"][b_idx]["b_fr"] = data["branch"][b_idx]["b_fr"] + parallel_br["b_fr"]
    data["branch"][b_idx]["b_to"] = data["branch"][b_idx]["b_to"] + parallel_br["b_to"]
    data["branch"][b_idx]["g_fr"] = data["branch"][b_idx]["g_fr"] + parallel_br["g_fr"]
    data["branch"][b_idx]["g_to"] = data["branch"][b_idx]["g_to"] + parallel_br["g_to"]
end


function generator_uc_data!(data; fcr_cost=50, droop_fac=1)
    for (g, gen) in data["gen"]
        gen["fcr_cost"] = fcr_cost * data["baseMVA"] / (3600 / data["frequency_parameters"]["t_fcrd"])  # 
        bus_id = gen["gen_bus"]
        if gen["gen_status"] == 1
            gen["zone"] = data["bus"]["$bus_id"]["zone"]
            gen["area"] = data["bus"]["$bus_id"]["area"]
            if haskey(gen, "H")
                gen["inertia_constants"] = parse(Float64, gen["H"])
            elseif gen["type"] == "Wind"
                gen["inertia_constants"] = 0.8
            else
                gen["inertia_constants"] = 0.0
            end
            if haskey(gen, "startup_warm_(dollar/MW)")
                gen["start_up_cost"] = gen["startup_warm_(dollar/MW)"] * data["baseMVA"]
            else
                gen["start_up_cost"] = 0.0
            end
            if haskey(gen, "Ramp_Up_Rate(MW/h)")
                gen["ramp_rate"] = (gen["Ramp_Up_Rate(MW/h)"] / data["baseMVA"]) / gen["pmax"] # in percent of maximum rating per hour
                gen["mdt"] = max(1, Int(gen["Minimum_Off_Time(Hours)"] / data["frequency_parameters"]["uc_time_interval"]))
                gen["mut"] = max(1, Int(gen["Minimum_On_Time(Hours)"] / data["frequency_parameters"]["uc_time_interval"]))
                gen["ramp_rate_per_s"] = (gen["Ramp_Up_Rate(MW/h)"] / data["baseMVA"] / 3600) * droop_fac
                gen["fcr_contribution"] = true
            else
                gen["ramp_rate"] = 1.0
                gen["ramp_rate_per_s"] = gen["pmax"] / 3600 * droop_fac
                gen["fcr_contribution"] = true
                gen["mdt"] = Int(1 / data["frequency_parameters"]["uc_time_interval"])
                gen["mut"] = Int(1 / data["frequency_parameters"]["uc_time_interval"])
            end
            idx = 1
            if haskey(gen, "variable_op_cost(dollar/MWh_sent_out)")
                if gen["type"] == "Fossil"
                    cost_var = max(25.0, gen["variable_op_cost(dollar/MWh_sent_out)"]) * data["baseMVA"]
                else
                    cost_var = gen["variable_op_cost(dollar/MWh_sent_out)"] * data["baseMVA"]
                end
                cost_fixed = ((gen["fixed_op_cost(dollar/MW/year)"]) * gen["pmax"] * data["baseMVA"] / (8760 / data["frequency_parameters"]["uc_time_interval"]))
                cost_fixed = 0
                gen["cost"] = [cost_var cost_fixed]
            elseif gen["type"] == "Wind" || gen["type"] == "Solar"
                cost = 5.0 + idx * 0.005
                gen["cost"] = [cost 0.0] * data["baseMVA"]
                idx += 1
            end
        else
            delete!(data["gen"], g)
        end
    end
end

function converter_uc_data!(data; t_hvdc=0.1, ffr_cost=50, rdc=1000)
    if haskey(data, "convdc")
        for (c, conv) in data["convdc"]
            conv_bus = conv["busac_i"]
            conv["zone"] = data["bus"]["$conv_bus"]["zone"]
            conv["area"] = data["bus"]["$conv_bus"]["area"]
            conv["t_hvdc"] = t_hvdc
            conv["ffr_cost"] = ffr_cost * data["baseMVA"] / (3600 / data["frequency_parameters"]["t_fcrd"])
            conv["rmax"] = rdc / data["baseMVA"]
        end
    end
end

function define_tie_lines!(data)
    data["tie_lines"] = Dict{String,Any}("1" => Dict(), "2" => Dict())
    # Heywood -> South East: 832 -> 1720
    # Dumeresq -> Bulli Creek:  104 -> 1623
    for (br, branch) in data["branch"]
        if (branch["f_bus"] == 832 && branch["t_bus"] == 1720)
            data["tie_lines"]["1"]["br_idx"] = parse(Int, br)
            data["tie_lines"]["1"]["area_fr"] = 1  # 1... NSW + VIC
            data["tie_lines"]["1"]["area_to"] = 4  # 4... SA 
        elseif (branch["t_bus"] == 832 && branch["f_bus"] == 1720)
            data["tie_lines"]["1"]["br_idx"] = parse(Int, br)
            data["tie_lines"]["1"]["area_to"] = 1  # 1... NSW + VIC
            data["tie_lines"]["1"]["area_fr"] = 4  # 4... SA 
        end
        if (branch["f_bus"] == 104 && branch["t_bus"] == 1623)
            data["tie_lines"]["2"]["br_idx"] = parse(Int, br)
            data["tie_lines"]["2"]["area_fr"] = 1  # 1... NSW + VIC
            data["tie_lines"]["2"]["area_to"] = 3  # 3... QLD 
        elseif (branch["t_bus"] == 104 && branch["f_bus"] == 1623)
            data["tie_lines"]["2"]["br_idx"] = parse(Int, br)
            data["tie_lines"]["2"]["area_to"] = 1  # 1... NSW + VIC
            data["tie_lines"]["2"]["area_fr"] = 3  # 3... QLD 
        end
    end
end







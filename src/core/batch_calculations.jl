function batch_fsuc(mn_data, fmin, dc_solver, scenario, year, h; droop = false, extension = "")

    objective_dc = Dict{String, Any}(["$i"=>[] for i in 1:length(fmin)])
    objective_no_dc = Dict{String, Any}(["$i"=>[] for i in 1:length(fmin)])
    t_dc = Dict{String, Any}(["$i"=>[] for i in 1:length(fmin)])
    t_no_dc = Dict{String, Any}(["$i"=>[] for i in 1:length(fmin)])
    t_opt_dc = Dict{String, Any}(["$i"=>[] for i in 1:length(fmin)])
    t_opt_no_dc = Dict{String, Any}(["$i"=>[] for i in 1:length(fmin)])

    y = "$year"
    if !isdir(joinpath("results", scenario))
        mkdir(joinpath("results", scenario))
    end
    if !isdir(joinpath("results", scenario, y))
        mkdir(joinpath("results", scenario, y))
    end
    if !isdir(joinpath("results", scenario, y, h))
        mkdir(joinpath("results", scenario, y, h))
    end
    
    filename = joinpath("results", scenario, y, h,  join(["input_data", extension,".json"]))
    json_string = JSON.json(mn_data)
    open(filename,"w") do f
    write(f, json_string)
    end
    json_string = nothing

    for idx = 1:length(fmin)
        fmin_ = fmin[idx]
        adjust_minimum_freqeuncy!(mn_data, fmin_)

        if droop == false
            s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "hvdc_inertia_contribution" => true, "relax_uc_binaries" => true)
            t_dc["$idx"] = @elapsed result_dc = CbaOPF.solve_fsuc(mn_data, _PM.DCPPowerModel, dc_solver, setting = s, multinetwork = true)
            filename = joinpath("results", scenario, y, h, join(["f",fmin_,extension,"_with_dc.json"]))
        else
            s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "hvdc_inertia_contribution" => true, "relax_uc_binaries" => true, "use_droop" => true)
            t_dc["$idx"] = @elapsed result_dc = CbaOPF.solve_fsuc(mn_data, _PM.DCPPowerModel, dc_solver, setting = s, multinetwork = true)
            filename = joinpath("results", scenario, y, h, join(["f",fmin_,extension,"_with_dc.json"]))
        end

        json_string = JSON.json(result_dc)
        open(filename,"w") do f
        write(f, json_string)
        end
        objective_dc["$idx"] = result_dc["objective"]
        t_opt_dc["$idx"] = result_dc["solve_time"]
        result_dc = nothing
        json_string = nothing

        if droop == false
            s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "hvdc_inertia_contribution" => false, "relax_uc_binaries" => true)
            t_no_dc["$idx"] = @elapsed result_no_dc = CbaOPF.solve_fsuc(mn_data, _PM.DCPPowerModel, dc_solver, setting = s, multinetwork = true)
            filename = joinpath("results", scenario, y, h, join(["f",fmin_,extension,"_without_dc.json"]))
        else
            s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "hvdc_inertia_contribution" => false, "relax_uc_binaries" => true, "use_droop" => true)
            t_no_dc["$idx"] = @elapsed result_no_dc = CbaOPF.solve_fsuc(mn_data, _PM.DCPPowerModel, dc_solver, setting = s, multinetwork = true)
            filename = joinpath("results", scenario, y, h, join(["f",fmin_,extension,"_without_dc.json"]))
        end
        json_string = JSON.json(result_no_dc)
        open(filename,"w") do f
        write(f, json_string)
        end
        objective_no_dc["$idx"] = result_no_dc["objective"]
        t_opt_no_dc["$idx"] = result_no_dc["solve_time"]
        result_no_dc = nothing
        json_string = nothing
    end

    filename = joinpath("results", scenario, y, h, join(["objective",extension,"_dc.json"]))
    json_string = JSON.json(objective_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    filename = joinpath("results", scenario, y, h, join(["objective",extension,"_no_dc.json"]))
    json_string = JSON.json(objective_no_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    filename = joinpath("results", scenario, y, h, join(["calculation",extension,"_time_dc.json"]))
    json_string = JSON.json(t_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    filename = joinpath("results", scenario, y, h, join(["calculation",extension,"_time_no_dc.json"]))
    json_string = JSON.json(t_no_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    filename = joinpath("results", scenario, y, h, join(["solver_time",extension,"_dc.json"]))
    json_string = JSON.json(t_opt_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    filename = joinpath("results", scenario, y, h, join(["solver_time",extension,"_no_dc.json"]))
    json_string = JSON.json(t_opt_no_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    return objective_dc, objective_no_dc, t_dc, t_no_dc
end

function adjust_minimum_freqeuncy!(data, fmin)
    if haskey(data, "nw")
        for (n, network) in data["nw"]
            network["frequency_parameters"]["fmin"] = fmin
            network["frequency_parameters"]["fmax"] =  network["frequency_parameters"]["f0"] + ((network["frequency_parameters"]["f0"] - fmin))
        end
    else
        data["frequency_parameters"]["fmin"] = fmin
        data["frequency_parameters"]["fmax"] = data["frequency_parameters"]["f0"] + ((data["frequency_parameters"]["f0"] - fmin))
    end
end



function batch_fsopf(data_dict, fmin, dc_solver, scenario, year, hours, h; droop = true, extension = extension)

    data = data_dict["data"]
    total_demand_series = data_dict["total_demand_series"] 
    dn_demand_series = data_dict["dn_demand_series"]
    pv_series = data_dict["pv_series"]
    wind_series = data_dict["wind_series"]
    wind_rez = data_dict["wind_rez"]
    pv_rez = data_dict["pv_rez"]
    generator_contingencies = data_dict["generator_contingencies"] 
    no_dc_cont = data_dict["no_dc_cont"]
    dn_res_factor = data_dict["dn_res_factor"]
    p2p = data_dict["p2p"] 


    objective_dc = Dict{String, Any}(["$i"=>[] for i in 1:length(fmin)])
    objective_no_dc = Dict{String, Any}(["$i"=>[] for i in 1:length(fmin)])
    t_dc = Dict{String, Any}(["$i"=>Dict{String, Any}(["$j"=>[] for j in 1:length(hours)]) for i in 1:length(fmin)])
    t_no_dc = Dict{String, Any}(["$i"=>Dict{String, Any}(["$j"=>[] for j in 1:length(hours)]) for i in 1:length(fmin)])
    t_opt_dc = Dict{String, Any}(["$i"=>Dict{String, Any}(["$j"=>[] for j in 1:length(hours)]) for i in 1:length(fmin)])
    t_opt_no_dc = Dict{String, Any}(["$i"=>Dict{String, Any}(["$j"=>[] for j in 1:length(hours)]) for i in 1:length(fmin)])
    
    y = "$year"
    if !isdir(joinpath("results", scenario))
        mkdir(joinpath("results", scenario))
    end
    if !isdir(joinpath("results", scenario, y))
        mkdir(joinpath("results", scenario, y))
    end
    if !isdir(joinpath("results", scenario, y, h))
        mkdir(joinpath("results", scenario, y, h))
    end
    
    for idx = 1:length(fmin)
        fmin_ = fmin[idx]
        adjust_minimum_freqeuncy!(data, fmin_)
    
        result_dc = Dict("objective" => 0.0, "solve_time" => 0.0, "solution" => Dict{String, Any}("nw" => Dict{String, Any}(["$i"=>Dict{String, Any}() for i in 1:length(hours)])))
        result_no_dc = Dict("objective" => 0.0, "solve_time" => 0.0, "solution" => Dict{String, Any}("nw" => Dict{String, Any}(["$i"=>Dict{String, Any}() for i in 1:length(hours)])))
        for h_idx in 1:length(hours)
            print("==================================== hour: ", h_idx, " ===============================================", "\n")
            print("=====================================================================================================", "\n")
            hour = hours[h_idx]
            mn_data = multi_network_uc_data(data, total_demand_series, dn_demand_series, pv_series, wind_series, pv_rez, wind_rez, hour, generator_contingencies, no_dc_cont = no_dc_cont, dn_res_factor = dn_res_factor, p2p = p2p)
            print("==================================== DC ====================================================", "\n")
            if droop == false
                s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "hvdc_inertia_contribution" => true, "relax_uc_binaries" => true)
                time_dc = @elapsed r_dc = CbaOPF.solve_fsuc(mn_data, _PM.DCPPowerModel, dc_solver, setting = s, multinetwork = true)
            else
                s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "hvdc_inertia_contribution" => true, "relax_uc_binaries" => true, "use_droop" => true)
                time_dc = @elapsed r_dc = CbaOPF.solve_fsuc(mn_data, _PM.DCPPowerModel, dc_solver, setting = s, multinetwork = true)
            end
            result_dc["objective"] = result_dc["objective"] + r_dc["objective"]
            result_dc["solution"]["nw"]["$h_idx"]  = r_dc["solution"]
    
            print("==================================== NO DC =================================================", "\n")
            if droop == false
                s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "hvdc_inertia_contribution" => false, "relax_uc_binaries" => true)
                time_no_dc = @elapsed r_no_dc = CbaOPF.solve_fsuc(mn_data, _PM.DCPPowerModel, dc_solver, setting = s, multinetwork = true)  
            else
                s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "hvdc_inertia_contribution" => false, "relax_uc_binaries" => true, "use_droop" => true)
                time_no_dc = @elapsed r_no_dc = CbaOPF.solve_fsuc(mn_data, _PM.DCPPowerModel, dc_solver, setting = s, multinetwork = true)
            end
            result_no_dc["objective"] = result_no_dc["objective"] + r_no_dc["objective"]
            result_no_dc["solution"]["nw"]["$h_idx"] = r_no_dc["solution"]
            t_dc["$idx"]["$h_idx"] = time_dc
            t_no_dc["$idx"]["$h_idx"] = time_no_dc
            t_opt_dc["$idx"]["$h_idx"] = r_dc["solve_time"]
            t_opt_no_dc["$idx"]["$h_idx"] = r_no_dc["solve_time"]
            mn_data = nothing
        end
    
        objective_dc["$idx"] = result_dc["objective"]
        objective_no_dc["$idx"] = result_no_dc["objective"]
    
        filename = joinpath("results", scenario, y, h, join(["f",fmin_,extension,"_with_dc.json"]))
        json_string = JSON.json(result_dc)
        open(filename,"w") do f
        write(f, json_string)
        end

        result_dc = nothing
        json_string = nothing
    
        filename = joinpath("results", scenario, y, h, join(["f",fmin_,extension,"_without_dc.json"]))
        json_string = JSON.json(result_no_dc)
        open(filename,"w") do f
        write(f, json_string)
        end
        result_no_dc = nothing
        json_string = nothing
    end
    
    filename = joinpath("results", scenario, y, h, join(["objective",extension,"_dc.json"]))
    json_string = JSON.json(objective_dc)
    open(filename,"w") do f
    write(f, json_string)
    end
    
    filename = joinpath("results", scenario, y, h, join(["objective",extension,"_no_dc.json"]))
    json_string = JSON.json(objective_no_dc)
    open(filename,"w") do f
    write(f, json_string)
    end
    
    filename = joinpath("results", scenario, y, h, join(["calculation",extension,"_time_dc.json"]))
    json_string = JSON.json(t_dc)
    open(filename,"w") do f
    write(f, json_string)
    end
    
    filename = joinpath("results", scenario, y, h, join(["calculation",extension,"_time_no_dc.json"]))
    json_string = JSON.json(t_no_dc)
    open(filename,"w") do f
    write(f, json_string)
    end
    
    filename = joinpath("results", scenario, y, h, join(["solver_time",extension,"_dc.json"]))
    json_string = JSON.json(t_opt_dc)
    open(filename,"w") do f
    write(f, json_string)
    end
    
    filename = joinpath("results", scenario, y, h, join(["solver_time",extension,"_no_dc.json"]))
    json_string = JSON.json(t_opt_no_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    return objective_dc, objective_no_dc, t_dc, t_no_dc, t_opt_dc, t_opt_no_dc
end
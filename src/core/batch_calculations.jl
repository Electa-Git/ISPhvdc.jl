function batch_fsuc(mn_data, fmin, dc_solver, scenario, year, h)

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

    filename = joinpath("results", scenario, y, h,  join(["input_data.json"]))
    json_string = JSON.json(mn_data)
    open(filename,"w") do f
    write(f, json_string)
    end
    json_string = nothing

    for idx = 1:length(fmin)
        fmin_ = fmin[idx]
        adjust_minimum_freqeuncy!(mn_data, fmin_)

        s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "hvdc_inertia_contribution" => true, "relax_uc_binaries" => true)
        t_dc["$idx"] = @elapsed result_dc = CbaOPF.solve_fsuc(mn_data, _PM.DCPPowerModel, dc_solver, setting = s, multinetwork = true)
        filename = joinpath("results", scenario, y, h, join(["f",fmin_,"_with_dc.json"]))
        json_string = JSON.json(result_dc)
        open(filename,"w") do f
        write(f, json_string)
        end
        objective_dc["$idx"] = result_dc["objective"]
        t_opt_dc["$idx"] = result_dc["solve_time"]
        result_dc = nothing
        json_string = nothing

        s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true, "hvdc_inertia_contribution" => false, "relax_uc_binaries" => true)
        t_no_dc["$idx"] = @elapsed result_no_dc = CbaOPF.solve_fsuc(mn_data, _PM.DCPPowerModel, dc_solver, setting = s, multinetwork = true)
        filename = joinpath("results", scenario, y, h, join(["f",fmin_,"_without_dc.json"]))
        json_string = JSON.json(result_no_dc)
        open(filename,"w") do f
        write(f, json_string)
        end
        objective_no_dc["$idx"] = result_no_dc["objective"]
        t_opt_no_dc["$idx"] = result_no_dc["solve_time"]
        result_no_dc = nothing
        json_string = nothing
    end

    filename = joinpath("results", scenario, y, h, join(["objective_dc.json"]))
    json_string = JSON.json(objective_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    filename = joinpath("results", scenario, y, h, join(["objective_no_dc.json"]))
    json_string = JSON.json(objective_no_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    filename = joinpath("results", scenario, y, h, join(["calculation_time_dc.json"]))
    json_string = JSON.json(t_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    filename = joinpath("results", scenario, y, h, join(["calculation_time_no_dc.json"]))
    json_string = JSON.json(t_no_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    filename = joinpath("results", scenario, y, h, join(["solver_time_dc.json"]))
    json_string = JSON.json(t_opt_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    filename = joinpath("results", scenario, y, h, join(["solver_time_no_dc.json"]))
    json_string = JSON.json(t_opt_no_dc)
    open(filename,"w") do f
    write(f, json_string)
    end

    return objective_dc, objective_no_dc, t_dc, t_no_dc
end

function adjust_minimum_freqeuncy!(mn_data, fmin)
    for (n, network) in mn_data["nw"]
        network["frequency_parameters"]["fmin"] = fmin
        network["frequency_parameters"]["fmax"] =  network["frequency_parameters"]["f0"] + ((network["frequency_parameters"]["f0"] - fmin))
    end
end
function plot_system_information(input_data, scenario, year, hours)
    hour_ids = input_data["hour_ids"]
    total_demand = zeros(length(hour_ids), 4)
    total_wind = zeros(length(hour_ids), 4)
    total_pv = zeros(length(hour_ids), 4)
    total_gen = zeros(length(hour_ids), 4)
    baseMVA = input_data["nw"]["1"]["baseMVA"]

    for idx in 1:length(hour_ids)
        hour = hour_ids[idx]

        for (l, load) in input_data["nw"]["$hour"]["load"]
            load_bus = load["load_bus"]
            area = input_data["nw"]["$hour"]["bus"]["$load_bus"]["area"]
            if area == 1
                total_demand[idx, 1] = total_demand[idx, 1]+ load["pd"] * baseMVA / 1e3
            elseif area == 3
                total_demand[idx, 2] = total_demand[idx, 2]+ load["pd"] * baseMVA / 1e3
            elseif area == 4
                total_demand[idx, 3] = total_demand[idx, 4]+ load["pd"] * baseMVA / 1e3
            elseif area == 5
                total_demand[idx, 4] = total_demand[idx, 4]+ load["pd"] * baseMVA / 1e3
            end
        end

        for (g, gen) in input_data["nw"]["$hour"]["gen"]
            gen_bus = gen["gen_bus"]
            area = input_data["nw"]["$hour"]["bus"]["$gen_bus"]["area"]
            if gen["type"] == "Wind"
                if area == 1
                    total_wind[idx, 1] = total_wind[idx, 1]+ gen["pmax"] * baseMVA / 1e3
                elseif area == 3
                    total_wind[idx, 2] = total_wind[idx, 2]+ gen["pmax"] * baseMVA / 1e3
                elseif area == 4
                    total_wind[idx, 3] = total_wind[idx, 4]+ gen["pmax"] * baseMVA / 1e3
                elseif area == 5
                    total_wind[idx, 4] = total_wind[idx, 4]+ gen["pmax"] * baseMVA / 1e3
                end
            elseif gen["type"] == "Solar"
                if area == 1
                    total_pv[idx, 1] = total_pv[idx, 1]+ gen["pmax"] * baseMVA / 1e3
                elseif area == 3
                    total_pv[idx, 2] = total_pv[idx, 2]+ gen["pmax"] * baseMVA / 1e3
                elseif area == 4
                    total_pv[idx, 3] = total_pv[idx, 4]+ gen["pmax"] * baseMVA / 1e3
                elseif area == 5
                    total_pv[idx, 4] = total_pv[idx, 4]+ gen["pmax"] * baseMVA / 1e3
                end
            end
            if area == 1
                total_gen[idx, 1] = total_gen[idx, 1]+ gen["pmax"] * gen["gen_status"] * baseMVA / 1e3
            elseif area == 3
                total_gen[idx, 2] = total_gen[idx, 2]+ gen["pmax"] * gen["gen_status"] * baseMVA / 1e3
            elseif area == 4
                total_gen[idx, 3] = total_gen[idx, 4]+ gen["pmax"] * gen["gen_status"] * baseMVA / 1e3
            elseif area == 5
                total_gen[idx, 4] = total_gen[idx, 4]+ gen["pmax"] * gen["gen_status"] * baseMVA / 1e3
            end
        end


    end

    legend = repeat(["NSW + VIC", "QLD", "SA", "TAS"], inner = length(hour_ids))
    x_values = repeat(["$(lpad(idx, 2, "0"))" for idx in 1:length(hour_ids)], outer = 4)
    p_d = StatsPlots.groupedbar(
        x_values, total_demand, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{d}~in~GW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["total_demand.pdf"]))
    StatsPlots.savefig(p_d, plot_filename)

    p_wind = StatsPlots.groupedbar(
        x_values, total_wind, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{w}~in~GW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["total_wind_capacity.pdf"]))
    StatsPlots.savefig(p_wind, plot_filename)

    p_pv = StatsPlots.groupedbar(
        x_values, total_pv, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{pv}~in~GW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["total_pv_capacity.pdf"]))
    StatsPlots.savefig(p_pv, plot_filename)

    p_gen = StatsPlots.groupedbar(
        x_values, total_gen, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{pv}~in~GW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    StatsPlots.plot!(p_gen, 0.5:1:length(hour_ids), sum(total_demand, dims = 2), marker = :diamond, label = "total demand", xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern")
    plot_filename = joinpath("results", scenario, year, hours,join(["total_gen_capacity.pdf"]))
    StatsPlots.savefig(p_gen, plot_filename)

end



function get_and_plot_objective_value(fmin, scenario, year, hours)

    fn = joinpath("results", scenario, year, hours, join(["objective_dc.json"]))
    objective_dc = Dict{String, Any}()
    open(fn) do f
        dicttxt = read(f,String)  # file information to string
        objective_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    fn = joinpath("results", scenario, year, hours, join(["objective_no_dc.json"]))
    objective_no_dc = Dict{String, Any}()
    open(fn) do f
        dicttxt = read(f,String)  # file information to string
        objective_no_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    o_dc = zeros(1, length(fmin))
    o_no_dc =  zeros(1, length(fmin))

    for idx in sort(parse.(Int, collect(keys(objective_dc))))
        if !isnothing(objective_dc["$idx"])
            o_dc[idx] = objective_dc["$idx"] / 1e6
        else
            o_dc[idx] = 0
        end

        if !isnothing(objective_no_dc["$idx"])
            o_no_dc[idx] = objective_no_dc["$idx"] / 1e6
        else
            o_no_dc[idx] = objective_no_dc["$(idx-1)"] / 1e6
        end
    end

    p1 = Plots.plot(fmin, o_no_dc', marker = :diamond, xlabel = "\$f_{min} in~Hz\$", ylabel = "\$Cost~in~M€\$", label = "without HVDC contribution", xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern")
    Plots.plot!(p1, fmin, o_dc', marker = :diamond,  xlabel = "\$f_{min} in~Hz\$", ylabel = "\$Cost~in~M€\$", label = "with HVDC contribution", xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern")
    plot_filename = joinpath("results", scenario, year, hours, "objective_comparison.pdf")
    Plots.savefig(p1, plot_filename)

    return objective_dc, objective_no_dc
end

function plot_calculation_time(fmin, scenario, year, hours)
    fn = joinpath("results", scenario, year, hours, join(["calculation_time_dc.json"]))
    time_dc = Dict{String, Any}()
    open(fn) do f
        dicttxt = read(f,String)  # file information to string
        time_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    fn = joinpath("results", scenario, year, hours, join(["calculation_time_no_dc.json"]))
    time_no_dc = Dict{String, Any}()
    open(fn) do f
        dicttxt = read(f,String)  # file information to string
        time_no_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    fn = joinpath("results", scenario, year, hours, join(["solver_time_dc.json"]))
    t_opt_dc = Dict{String, Any}()
    open(fn) do f
        dicttxt = read(f,String)  # file information to string
        t_opt_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    fn = joinpath("results", scenario, year, hours, join(["solver_time_no_dc.json"]))
    t_opt_no_dc = Dict{String, Any}()
    open(fn) do f
        dicttxt = read(f,String)  # file information to string
        t_opt_no_dc = JSON.parse(dicttxt)  # parse and transform data
    end


    t_no_dc = zeros(length(fmin), 2)
    t_dc = zeros(length(fmin), 2)

    for idx in 1:length(fmin)
        t_no_dc[idx, 1] = t_opt_no_dc["$idx"]
        t_no_dc[idx, 2] = time_no_dc["$idx"] - t_opt_no_dc["$idx"]

        t_dc[idx, 1] = t_opt_dc["$idx"]
        t_dc[idx, 2] = time_dc["$idx"] - t_opt_dc["$idx"]
    end

    legend = repeat(["Solver time", "Model building"], inner = length(fmin))
    x_values = repeat(["$f" for f in fmin], outer = 2)

    p_no_dc = StatsPlots.groupedbar(
        x_values, t_no_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$f_{min}~in~Hz\$", ylabel = "\$t~in~s\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["calculation_time_no_dc.pdf"]))
    StatsPlots.savefig(p_no_dc, plot_filename)

    p_dc = StatsPlots.groupedbar(
        x_values, t_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$f_{min}~in~Hz\$", ylabel = "\$t~in~s\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["calculation_time_dc.pdf"]))
    StatsPlots.savefig(p_dc, plot_filename)
end

function load_input_data(scenario, year, hours)
    print("Loading input data", "\n")
    
    fn = joinpath("results", scenario, year, hours, join(["input_data.json"]))
    input_data = Dict{String, Any}()
    open(fn) do f
        dicttxt = read(f,String)  # file information to string
        input_data = JSON.parse(dicttxt)  # parse and transform data
    end

    return input_data
end

function plot_load_shedding(input_data, fmin, scenario, year, hours)
    print("Loading results", "\n")
    fn = joinpath("results", scenario, year, hours, join(["f",fmin,"_with_dc.json"]))
    result_dc = Dict{String, Any}()
    open(fn) do f
    dicttxt = read(f,String)  # file information to string
        result_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    fn = joinpath("results", scenario, year, hours, join(["f",fmin,"_without_dc.json"]))
    result_no_dc = Dict{String, Any}()
    open(fn) do f
    dicttxt = read(f,String)  # file information to string
        result_no_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    hour_ids = input_data["hour_ids"]

    ls_no_dc = zeros(length(hour_ids), 4)
    ls_dc = zeros(length(hour_ids), 4)

    baseMVA = input_data["nw"]["1"]["baseMVA"]

    print("Determining demand shedding", "\n")
    for idx in 1:length(hour_ids)
        hour = hour_ids[idx]

        res_no_dc = result_no_dc["solution"]["nw"]["$hour"]
        for (l, load) in res_no_dc["load"]
            load_bus = input_data["nw"]["1"]["load"][l]["load_bus"]
            area = input_data["nw"]["1"]["bus"]["$load_bus"]["area"]
            if area == 1
                ls_no_dc[idx, 1] = ls_no_dc[idx, 1]+ load["pcurt"] * baseMVA
            elseif area == 3
                ls_no_dc[idx, 2] = ls_no_dc[idx, 2] + load["pcurt"] * baseMVA
            elseif area == 4
                ls_no_dc[idx, 3] = ls_no_dc[idx, 3] + load["pcurt"] * baseMVA
            elseif area == 5
                ls_no_dc[idx, 4] = ls_no_dc[idx, 4] + load["pcurt"] * baseMVA
            end
        end

        res_dc = result_dc["solution"]["nw"]["$hour"]
        for (l, load) in res_dc["load"]
            load_bus = input_data["nw"]["1"]["load"][l]["load_bus"]
            area = input_data["nw"]["1"]["bus"]["$load_bus"]["area"]
            if area == 1
                ls_dc[idx, 1] = ls_dc[idx, 1]+ load["pcurt"] * baseMVA
            elseif area == 3
                ls_dc[idx, 2] = ls_dc[idx, 2] + load["pcurt"] * baseMVA
            elseif area == 4
                ls_dc[idx, 3] = ls_dc[idx, 3] + load["pcurt"] * baseMVA
            elseif area == 5
                ls_dc[idx, 4] = ls_dc[idx, 4] + load["pcurt"] * baseMVA
            end
        end
    end

    print("Plotting", "\n")

    legend = repeat(["NSW + VIC w/o dc", "QLD w/o dc", "SA w/o dc", "TAS w/o dc"], inner = length(hour_ids))
    x_values = repeat(["$(lpad(idx, 2, "0"))" for idx in 1:length(hour_ids)], outer = 4)
    p_no_dc = StatsPlots.groupedbar(
        x_values, ls_no_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P^{curt}~in~MW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours, join(["demand_shedding_without_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_no_dc, plot_filename)

    legend = repeat(["NSW + VIC with dc", "QLD with dc", "SA with dc", "TAS with dc"], inner = length(hour_ids))
    p_dc = StatsPlots.groupedbar(
        x_values, ls_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P^{curt}~in~MW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours, join(["demand_shedding_with_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_dc, plot_filename)






    # p_comp = Plots.plot(1:length(hour_ids), ls_no_dc[:, 1], linestyle = :solid, linecolor = :black, label = "NSW + VIC, no dc", marker = :diamond, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright)
    # Plots.plot!(p_comp, 1:length(hour_ids), ls_dc[:, 1], linestyle = :dash, linecolor = :black, label = "NSW + VIC, with dc", marker = :diamond, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright)
    # Plots.plot!(p_comp, 1:length(hour_ids), ls_no_dc[:, 2], linestyle = :solid, linecolor = :black, label = "QLD, no dc", marker = :circle, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright)
    # Plots.plot!(p_comp, 1:length(hour_ids), ls_dc[:, 2], linestyle = :dash, linecolor = :black, label = "QLD, with dc", marker= :circle, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright)
    # Plots.plot!(p_comp, 1:length(hour_ids), ls_no_dc[:, 3], linestyle = :solid, linecolor = :black, label = "SA, no dc", marker = :square, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright)
    # Plots.plot!(p_comp, 1:length(hour_ids), ls_dc[:, 3], linestyle = :dash, linecolor = :black, label = "SA, with dc", marker = :square, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright)
    # Plots.plot!(p_comp, 1:length(hour_ids), ls_no_dc[:, 4], linestyle = :solid, linecolor = :black, label = "TAS, no dc", marker = :xcross, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright)
    # Plots.plot!(p_comp, 1:length(hour_ids), ls_dc[:, 4], linestyle = :dash, linecolor = :black, label = "TAS, with dc", marker = :xcross, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright)

    # plot_filename = joinpath("results", scenario, year, hours, join(["demandshedding_f",fmin,".pdf"]))
    # StatsPlots.savefig(p_comp, plot_filename)



end

function plot_res_generation_and_curtailment(input_data, fmin, scenario, year, hours)
    print("Loading results", "\n")
    fn = joinpath("results", scenario, year, hours, join(["f",fmin,"_with_dc.json"]))
    result_dc = Dict{String, Any}()
    open(fn) do f
    dicttxt = read(f,String)  # file information to string
        result_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    fn = joinpath("results", scenario, year, hours, join(["f",fmin,"_without_dc.json"]))
    result_no_dc = Dict{String, Any}()
    open(fn) do f
    dicttxt = read(f,String)  # file information to string
        result_no_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    hour_ids = input_data["hour_ids"]

    curt_pv_no_dc = zeros(length(hour_ids), 4)
    curt_pv_dc = zeros(length(hour_ids), 4)
    gen_pv_dc = zeros(length(hour_ids), 4)
    gen_pv_no_dc = zeros(length(hour_ids), 4)

    curt_wind_no_dc = zeros(length(hour_ids), 4)
    curt_wind_dc = zeros(length(hour_ids), 4)
    gen_wind_dc = zeros(length(hour_ids), 4)
    gen_wind_no_dc = zeros(length(hour_ids), 4)

    baseMVA = input_data["nw"]["1"]["baseMVA"]

    print("Determining RES curtailment", "\n")
    for idx in 1:length(hour_ids)
        hour = hour_ids[idx]

        res_no_dc = result_no_dc["solution"]["nw"]["$hour"]
        for (g, gen) in res_no_dc["gen"]
            gen_bus = input_data["nw"]["1"]["gen"][g]["gen_bus"]
            area = input_data["nw"]["1"]["bus"]["$gen_bus"]["area"]
            if input_data["nw"]["$hour"]["gen"][g]["type"] == "Wind"
                pmax = input_data["nw"]["$hour"]["gen"][g]["pmax"]
                if area == 1
                    curt_wind_no_dc[idx, 1] = curt_wind_no_dc[idx, 1]+ (pmax - gen["pg"]) * baseMVA
                    gen_wind_no_dc[idx, 1] = gen_wind_no_dc[idx, 1]+ gen["pg"] * baseMVA
                elseif area == 3
                    curt_wind_no_dc[idx, 2] = curt_wind_no_dc[idx, 2]+ (pmax - gen["pg"]) * baseMVA
                    gen_wind_no_dc[idx, 2] = gen_wind_no_dc[idx, 2]+ gen["pg"] * baseMVA
                elseif area == 4
                    curt_wind_no_dc[idx, 3] = curt_wind_no_dc[idx, 3]+ (pmax - gen["pg"]) * baseMVA
                    gen_wind_no_dc[idx, 3] = gen_wind_no_dc[idx, 3]+ gen["pg"] * baseMVA
                elseif area == 5
                    curt_wind_no_dc[idx, 4] = curt_wind_no_dc[idx, 4]+ (pmax - gen["pg"]) * baseMVA
                    gen_wind_no_dc[idx, 4] = gen_wind_no_dc[idx, 4]+ gen["pg"] * baseMVA
                end
            elseif input_data["nw"]["$hour"]["gen"][g]["type"] == "Solar"
                pmax = input_data["nw"]["$hour"]["gen"][g]["pmax"]
                if area == 1
                    curt_pv_no_dc[idx, 1] = curt_pv_no_dc[idx, 1]+ (pmax - gen["pg"]) * baseMVA
                    gen_pv_no_dc[idx, 1] = gen_pv_no_dc[idx, 1]+ gen["pg"] * baseMVA
                elseif area == 3
                    curt_pv_no_dc[idx, 2] = curt_pv_no_dc[idx, 2]+ (pmax - gen["pg"]) * baseMVA
                    gen_pv_no_dc[idx, 2] = gen_pv_no_dc[idx, 2]+ gen["pg"] * baseMVA
                elseif area == 4
                    curt_pv_no_dc[idx, 3] = curt_pv_no_dc[idx, 3]+ (pmax - gen["pg"]) * baseMVA
                    gen_pv_no_dc[idx, 3] = gen_pv_no_dc[idx, 3]+ gen["pg"] * baseMVA
                elseif area == 5
                    curt_pv_no_dc[idx, 4] = curt_pv_no_dc[idx, 4]+ (pmax - gen["pg"]) * baseMVA
                    gen_pv_no_dc[idx, 4] = gen_pv_no_dc[idx, 4]+ gen["pg"] * baseMVA
                end
            end
        end

        res_dc = result_dc["solution"]["nw"]["$hour"]
        for (g, gen) in res_dc["gen"]
            gen_bus = input_data["nw"]["1"]["gen"][g]["gen_bus"]
            area = input_data["nw"]["1"]["bus"]["$gen_bus"]["area"]
            if input_data["nw"]["$hour"]["gen"][g]["type"] == "Wind"
                pmax = input_data["nw"]["$hour"]["gen"][g]["pmax"]
                if area == 1
                    curt_wind_dc[idx, 1] = curt_wind_dc[idx, 1]+ (pmax - gen["pg"]) * baseMVA
                    gen_wind_dc[idx, 1] = gen_wind_dc[idx, 1]+ gen["pg"] * baseMVA
                elseif area == 3
                    curt_wind_dc[idx, 2] = curt_wind_dc[idx, 2]+ (pmax - gen["pg"]) * baseMVA
                    gen_wind_dc[idx, 2] = gen_wind_dc[idx, 2]+ gen["pg"] * baseMVA
                elseif area == 4
                    curt_wind_dc[idx, 3] = curt_wind_dc[idx, 3]+ (pmax - gen["pg"]) * baseMVA
                    gen_wind_dc[idx, 3] = gen_wind_dc[idx, 3]+ gen["pg"] * baseMVA
                elseif area == 5
                    curt_wind_dc[idx, 4] = curt_wind_dc[idx, 4]+ (pmax - gen["pg"]) * baseMVA
                    gen_wind_dc[idx, 4] = gen_wind_dc[idx, 4]+ gen["pg"] * baseMVA
                end
            elseif input_data["nw"]["$hour"]["gen"][g]["type"] == "Solar"
                pmax = input_data["nw"]["$hour"]["gen"][g]["pmax"]
                if area == 1
                    curt_pv_dc[idx, 1] = curt_pv_dc[idx, 1]+ (pmax - gen["pg"]) * baseMVA
                    gen_pv_dc[idx, 1] = gen_pv_dc[idx, 1]+ gen["pg"] * baseMVA
                elseif area == 3
                    curt_pv_dc[idx, 2] = curt_pv_dc[idx, 2]+ (pmax - gen["pg"]) * baseMVA
                    gen_pv_dc[idx, 2] = gen_pv_dc[idx, 2]+ gen["pg"] * baseMVA
                elseif area == 4
                    curt_pv_dc[idx, 3] = curt_pv_dc[idx, 3]+ (pmax - gen["pg"]) * baseMVA
                    gen_pv_dc[idx, 3] = gen_pv_dc[idx, 3]+ gen["pg"] * baseMVA
                elseif area == 5
                    curt_pv_dc[idx, 4] = curt_pv_dc[idx, 4]+ (pmax - gen["pg"]) * baseMVA
                    gen_pv_dc[idx, 4] = gen_pv_dc[idx, 4]+ gen["pg"] * baseMVA
                end
            end
        end
        end

    print("Plotting", "\n")

    legend = repeat(["NSW + VIC", "QLD", "SA", "TAS"], inner = length(hour_ids))
    x_values = repeat(["$(lpad(idx, 2, "0"))" for idx in 1:length(hour_ids)], outer = 4)
    p_wind_no_dc = StatsPlots.groupedbar(
        x_values, curt_wind_no_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{curt}~in~MW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["wind_curtailment_without_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_wind_no_dc, plot_filename)

    p_wind_dc = StatsPlots.groupedbar(
        x_values, curt_wind_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{curt}~in~MW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["wind_curtailment_with_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_wind_dc, plot_filename)

    p_pv_no_dc = StatsPlots.groupedbar(
        x_values, curt_pv_no_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{curt}~in~MW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["pv_curtailment_without_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_pv_no_dc, plot_filename)

    p_pv_dc = StatsPlots.groupedbar(
        x_values, curt_pv_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{curt}~in~MW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["pv_curtailment_with_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_pv_dc, plot_filename)

    p_wind_no_dc = StatsPlots.groupedbar(
        x_values, gen_wind_no_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{g}~in~MW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["wind_generation_without_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_wind_no_dc, plot_filename)

    p_wind_dc = StatsPlots.groupedbar(
        x_values, gen_wind_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{g}~in~MW\$", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["wind_generation_with_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_wind_dc, plot_filename)

    p_pv_no_dc = StatsPlots.groupedbar(
        x_values, gen_pv_no_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{g}~in~MW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["pv_generation_without_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_pv_no_dc, plot_filename)

    p_pv_dc = StatsPlots.groupedbar(
        x_values, gen_pv_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$P_{g}~in~MW\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours,join(["pv_generation_with_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_pv_dc, plot_filename)
end


function plot_total_inertia(input_data, fmin, scenario, year, hours)
    print("Loading results", "\n")
    fn = joinpath("results", scenario, year, hours, join(["f",fmin,"_with_dc.json"]))
    result_dc = Dict{String, Any}()
    open(fn) do f
    dicttxt = read(f,String)  # file information to string
        result_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    fn = joinpath("results", scenario, year, hours, join(["f",fmin,"_without_dc.json"]))
    result_no_dc = Dict{String, Any}()
    open(fn) do f
    dicttxt = read(f,String)  # file information to string
        result_no_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    hour_ids = input_data["hour_ids"]

    in_no_dc = zeros(length(hour_ids), 4)
    in_dc = zeros(length(hour_ids), 4)

    baseMVA = input_data["nw"]["1"]["baseMVA"]

    print("Determining inertia", "\n")
    for idx in 1:length(hour_ids)
        hour = hour_ids[idx]

        res_no_dc = result_no_dc["solution"]["nw"]["$hour"]
        for (g, gen) in res_no_dc["gen"]
            gen_bus = input_data["nw"]["1"]["gen"][g]["gen_bus"]
            area = input_data["nw"]["1"]["bus"]["$gen_bus"]["area"]
            if area == 1
                in_no_dc[idx, 1] = in_no_dc[idx, 1]+ input_data["nw"]["$hour"]["gen"][g]["pmax"] * input_data["nw"]["$hour"]["gen"][g]["inertia_constants"] * gen["alpha_g"] * baseMVA / 1e3
            elseif area == 3
                in_no_dc[idx, 2] = in_no_dc[idx, 2]+ input_data["nw"]["$hour"]["gen"][g]["pmax"] * input_data["nw"]["$hour"]["gen"][g]["inertia_constants"] * gen["alpha_g"] * baseMVA / 1e3
            elseif area == 4
                in_no_dc[idx, 3] = in_no_dc[idx, 3]+ input_data["nw"]["$hour"]["gen"][g]["pmax"] * input_data["nw"]["$hour"]["gen"][g]["inertia_constants"] * gen["alpha_g"] * baseMVA / 1e3
            elseif area == 5
                in_no_dc[idx, 4] = in_no_dc[idx, 4]+ input_data["nw"]["$hour"]["gen"][g]["pmax"] * input_data["nw"]["$hour"]["gen"][g]["inertia_constants"] * gen["alpha_g"] * baseMVA / 1e3
            end
        end

        res_dc = result_dc["solution"]["nw"]["$hour"]
        for (g, gen) in res_dc["gen"]
            gen_bus = input_data["nw"]["1"]["gen"][g]["gen_bus"]
            area = input_data["nw"]["1"]["bus"]["$gen_bus"]["area"]
            if area == 1
                in_dc[idx, 1] = in_dc[idx, 1]+ input_data["nw"]["$hour"]["gen"][g]["pmax"] * input_data["nw"]["$hour"]["gen"][g]["inertia_constants"] * gen["alpha_g"] * baseMVA / 1e3
            elseif area == 3
                in_dc[idx, 2] = in_dc[idx, 2]+ input_data["nw"]["$hour"]["gen"][g]["pmax"] * input_data["nw"]["$hour"]["gen"][g]["inertia_constants"] * gen["alpha_g"] * baseMVA / 1e3
            elseif area == 4
                in_dc[idx, 3] = in_dc[idx, 3]+ input_data["nw"]["$hour"]["gen"][g]["pmax"] * input_data["nw"]["$hour"]["gen"][g]["inertia_constants"] * gen["alpha_g"] * baseMVA / 1e3
            elseif area == 5
                in_dc[idx, 4] = in_dc[idx, 4]+ input_data["nw"]["$hour"]["gen"][g]["pmax"] * input_data["nw"]["$hour"]["gen"][g]["inertia_constants"] * gen["alpha_g"] * baseMVA / 1e3
            end
        end
    end

    print("Plotting", "\n")

    legend = repeat(["NSW + VIC", "QLD", "SA", "TAS"], inner = length(hour_ids))
    x_values = repeat(["$(lpad(idx, 2, "0"))" for idx in 1:length(hour_ids)], outer = 4)
    p_no_dc = StatsPlots.groupedbar(
        x_values, in_no_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours, join(["inertia_without_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_no_dc, plot_filename)

    p_dc = StatsPlots.groupedbar(
        x_values, in_dc, group = legend,
        bar_position = :stack,
        xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$",
        xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern"
    )
    plot_filename = joinpath("results", scenario, year, hours, join(["inertia_with_dc_f",fmin,".pdf"]))
    StatsPlots.savefig(p_dc, plot_filename)


    p_comp = Plots.plot(1:length(hour_ids), in_no_dc[:, 1], linestyle = :solid, linecolor = :black, label = "NSW + VIC, no dc", marker = :diamond, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright, fontfamily = "Computer Modern")
    Plots.plot!(p_comp, 1:length(hour_ids), in_dc[:, 1], linestyle = :dash, linecolor = :black, label = "NSW + VIC, with dc", marker = :diamond, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright, fontfamily = "Computer Modern")
    Plots.plot!(p_comp, 1:length(hour_ids), in_no_dc[:, 2], linestyle = :solid, linecolor = :black, label = "QLD, no dc", marker = :circle, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright, fontfamily = "Computer Modern")
    Plots.plot!(p_comp, 1:length(hour_ids), in_dc[:, 2], linestyle = :dash, linecolor = :black, label = "QLD, with dc", marker= :circle, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright, fontfamily = "Computer Modern")
    Plots.plot!(p_comp, 1:length(hour_ids), in_no_dc[:, 3], linestyle = :solid, linecolor = :black, label = "SA, no dc", marker = :square, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright, fontfamily = "Computer Modern")
    Plots.plot!(p_comp, 1:length(hour_ids), in_dc[:, 3], linestyle = :dash, linecolor = :black, label = "SA, with dc", marker = :square, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright, fontfamily = "Computer Modern")
    Plots.plot!(p_comp, 1:length(hour_ids), in_no_dc[:, 4], linestyle = :solid, linecolor = :black, label = "TAS, no dc", marker = :xcross, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright, fontfamily = "Computer Modern")
    Plots.plot!(p_comp, 1:length(hour_ids), in_dc[:, 4], linestyle = :dash, linecolor = :black, label = "TAS, with dc", marker = :xcross, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$H~in~GWs\$", legend_position = :outertopright, fontfamily = "Computer Modern")

    plot_filename = joinpath("results", scenario, year, hours, join(["inertia_comparison_f",fmin,".pdf"]))
    StatsPlots.savefig(p_comp, plot_filename)
end


function plot_tie_line_flows(input_data, fmin, scenario, year, hours)
    print("Loading results", "\n")
    fn = joinpath("results", scenario, year, hours, join(["f",fmin,"_with_dc.json"]))
    result_dc = Dict{String, Any}()
    open(fn) do f
    dicttxt = read(f,String)  # file information to string
        result_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    fn = joinpath("results", scenario, year, hours, join(["f",fmin,"_without_dc.json"]))
    result_no_dc = Dict{String, Any}()
    open(fn) do f
    dicttxt = read(f,String)  # file information to string
        result_no_dc = JSON.parse(dicttxt)  # parse and transform data
    end

    hour_ids = input_data["hour_ids"]

    tlf_no_dc = zeros(length(hour_ids), length(input_data["nw"]["1"]["tie_lines"]))
    tlf_dc = zeros(length(hour_ids), length(input_data["nw"]["1"]["tie_lines"]))

    baseMVA = input_data["nw"]["1"]["baseMVA"]

    print("Determining tie line flows", "\n")
    for idx in 1:length(hour_ids)
        hour = hour_ids[idx]

        res_no_dc = result_no_dc["solution"]["nw"]["$hour"]
        for (b, branch) in res_no_dc["branch"]
            for (t, tl) in input_data["nw"]["1"]["tie_lines"]
                if parse(Int, b) == tl["br_idx"]
                    tlf_no_dc[idx, parse(Int, t)] = branch["pf"] * baseMVA
                end
            end
        end

        res_dc = result_dc["solution"]["nw"]["$hour"]
        for (b, branch) in res_dc["branch"]
            for (t, tl) in input_data["nw"]["1"]["tie_lines"]
                if parse(Int, b) == tl["br_idx"]
                    tlf_dc[idx, parse(Int, t)] = branch["pf"] * baseMVA
                end
            end
        end
    end


    print("Plotting", "\n")

    for (t, tl) in input_data["nw"]["1"]["tie_lines"]
        if tl["area_fr"] == 1
            area_fr = "NSW + VIC"
        elseif tl["area_fr"] == 3
            area_fr = "QLD"
        elseif tl["area_fr"] == 4
            area_fr = "SA"
        elseif tl["area_fr"] == 5
            area_fr = "TAS"
        end
        if tl["area_to"] == 1
            area_to = "NSW + VIC"
        elseif tl["area_to"] == 3
            area_to = "QLD"
        elseif tl["area_to"] == 4
            area_to = "SA"
        elseif tl["area_to"] == 5
            area_to = "TAS"
        end
        tl["direction"] = join([area_fr, " -> ", area_to])
    end



    label = join([input_data["nw"]["1"]["tie_lines"]["1"]["direction"], " no dc"])    
    p_comp = Plots.plot(1:length(hour_ids), tlf_no_dc[:, 1], linestyle = :solid, label = label, marker = :diamond, xlabel = "\$hour~id\$", ylabel = "\$P_{tl}~in~MVA\$", legend_position = :outertopright,xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern")
    label = join([input_data["nw"]["1"]["tie_lines"]["1"]["direction"], " with dc"])  
    Plots.plot!(p_comp, 1:length(hour_ids), tlf_dc[:, 1], linestyle = :solid, label = label, marker = :diamond, xlabel = "\$hour~id\$", ylabel = "\$P_{tl}~in~MVA\$", legend_position = :outertopright,xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern")
    for idx in 2:length(input_data["nw"]["1"]["tie_lines"])
        label = join([input_data["nw"]["1"]["tie_lines"]["$idx"]["direction"], " no dc"])    
        Plots.plot!(p_comp, 1:length(hour_ids), tlf_no_dc[:, idx], linestyle = :solid, label = label, marker = :diamond, xlabel = "\$hour~id\$", ylabel = "\$P_{tl}~in~MVA\$", legend_position = :outertopright,xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern")
        label = join([input_data["nw"]["1"]["tie_lines"]["$idx"]["direction"], " with dc"])  
        Plots.plot!(p_comp, 1:length(hour_ids), tlf_dc[:, idx], linestyle = :solid, label = label, marker = :diamond, xlabel = "\$hour~id\$", ylabel = "\$P_{tl}~in~MVA\$", legend_position = :outertopright,xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern")
    end

    plot_filename = joinpath("results", scenario, year, hours, join(["flow_comparison_f",fmin,".pdf"]))
    StatsPlots.savefig(p_comp, plot_filename)
end


# function plot_tie_line_flows(input_data, fmin, scenario, year, hours)
#     print("Loading results", "\n")
#     fn = joinpath("results", scenario, year, hours, join(["f",fmin,"_with_dc.json"]))
#     result_dc = Dict{String, Any}()
#     open(fn) do f
#     dicttxt = read(f,String)  # file information to string
#         result_dc = JSON.parse(dicttxt)  # parse and transform data
#     end

#     fn = joinpath("results", scenario, year, hours, join(["f",fmin,"_without_dc.json"]))
#     result_no_dc = Dict{String, Any}()
#     open(fn) do f
#     dicttxt = read(f,String)  # file information to string
#         result_no_dc = JSON.parse(dicttxt)  # parse and transform data
#     end

#     hour_ids = input_data["hour_ids"]

#     dcf_no_dc = zeros(length(hour_ids), length(input_data["nw"]["1"]["branchdc"]))
#     dcf_dc = zeros(length(hour_ids), length(input_data["nw"]["1"]["branchdc"]))

#     baseMVA = input_data["nw"]["1"]["baseMVA"]

#     print("Determining tie line flows", "\n")
#     for idx in 1:length(hour_ids)
#         hour = hour_ids[idx]

#         res_no_dc = result_no_dc["solution"]["nw"]["$hour"]
#         for (b, branch) in res_no_dc["branchdc"]
#             dcf_no_dc[idx, parse(Int, t)] = branch["pf"] * baseMVA
#         end

#         res_dc = result_dc["solution"]["nw"]["$hour"]
#         for (b, branch) in res_dc["branch"]
#             dcf_dc[idx, parse(Int, t)] = branch["pf"] * baseMVA
#         end
#     end


#     print("Plotting", "\n")

#     for (t, tl) in input_data["nw"]["1"]["tie_lines"]
#         if tl["area_fr"] == 1
#             area_fr = "NSW + VIC"
#         elseif tl["area_fr"] == 3
#             area_fr = "QLD"
#         elseif tl["area_fr"] == 4
#             area_fr = "SA"
#         elseif tl["area_fr"] == 5
#             area_fr = "TAS"
#         end
#         if tl["area_to"] == 1
#             area_to = "NSW + VIC"
#         elseif tl["area_to"] == 3
#             area_to = "QLD"
#         elseif tl["area_to"] == 4
#             area_to = "SA"
#         elseif tl["area_to"] == 5
#             area_to = "TAS"
#         end
#         tl["direction"] = join([area_fr, " -> ", area_to])
#     end



#     label = join([input_data["nw"]["1"]["tie_lines"]["1"]["direction"], " no dc"])    
#     p_comp = Plots.plot(1:length(hour_ids), tlf_no_dc[:, 1], linestyle = :solid, label = label, marker = :diamond, xlabel = "\$hour~id\$", ylabel = "\$P_{tl}~in~MVA\$", legend_position = :outertopright)
#     label = join([input_data["nw"]["1"]["tie_lines"]["1"]["direction"], " with dc"])  
#     Plots.plot!(p_comp, 1:length(hour_ids), tlf_dc[:, 1], linestyle = :solid, label = label, marker = :diamond, xlabel = "\$hour~id\$", ylabel = "\$P_{tl}~in~MVA\$", legend_position = :outertopright)
#     for idx in 2:length(input_data["nw"]["1"]["tie_lines"])
#         label = join([input_data["nw"]["1"]["tie_lines"]["$idx"]["direction"], " no dc"])    
#         Plots.plot!(p_comp, 1:length(hour_ids), tlf_no_dc[:, idx], linestyle = :solid, label = label, marker = :diamond, xlabel = "\$hour~id\$", ylabel = "\$P_{tl}~in~MVA\$", legend_position = :outertopright)
#         label = join([input_data["nw"]["1"]["tie_lines"]["$idx"]["direction"], " with dc"])  
#         Plots.plot!(p_comp, 1:length(hour_ids), tlf_dc[:, idx], linestyle = :solid, label = label, marker = :diamond, xlabel = "\$hour~id\$", ylabel = "\$P_{tl}~in~MVA\$", legend_position = :outertopright)
#     end

#     plot_filename = joinpath("results", scenario, year, hours, join(["flow_comparison_f",fmin,".pdf"]))
#     StatsPlots.savefig(p_comp, plot_filename)
# end


function plot_hvdc_contribution(input_data, fmin_, scenario, year, hours)
    p_in_hvdc = Plots.plot()
    for fmin in fmin_
        print("Loading results", "\n")
        fn = joinpath("results", scenario, year, hours, join(["f",fmin,"_with_dc.json"]))
        result_dc = Dict{String, Any}()
        open(fn) do f
        dicttxt = read(f,String)  # file information to string
            result_dc = JSON.parse(dicttxt)  # parse and transform data
        end

        hour_ids = input_data["hour_ids"]

        in_dc = zeros(length(hour_ids), input_data["number_of_contingencies"])

        baseMVA = input_data["nw"]["1"]["baseMVA"]

        print("Determining inertia", "\n")
        for idx in 1:length(hour_ids)
            hour = hour_ids[idx]
            conts = hour .+ collect(1:(input_data["number_of_contingencies"] -1))
            for jdx in 1:length(conts)
                cont = conts[jdx]
                in_dc[idx, jdx] = sum([conv["pconv_in_abs"] for (c, conv)  in  result_dc["solution"]["nw"]["$cont"]["convdc"]]) * baseMVA / 1e3
                if in_dc[idx, jdx] !== 0.0
                    cont_name = input_data["nw"]["$cont"]["contingency"]
                    if !isnothing(cont_name["gen_id"])
                        gen_id = cont_name["gen_id"]
                        cont_mw = result_dc["solution"]["nw"]["$hour"]["gen"]["$gen_id"]["pg"] * input_data["nw"]["1"]["baseMVA"]
                        cont_string = join(["generator ", gen_id, " with ", cont_mw, " MW"])
                    elseif !isnothing(cont_name["conv_id"])
                        conv_id = cont_name["conv_id"]
                        cont_mw = result_dc["solution"]["nw"]["$hour"]["convdc"]["$conv_id"]["pgrid"] * input_data["nw"]["1"]["baseMVA"]
                        cont_string = join(["converter ", conv_id, " with ", cont_mw, " MW"])
                    elseif !isnothing(cont_name["branch_id"])
                        branch_id = input_data["nw"]["$hour"]["tie_lines"]["$(cont_name["branch_id"])"]["br_idx"]
                        cont_mw = result_dc["solution"]["nw"]["$hour"]["branch"]["$branch_id"]["pf"] * input_data["nw"]["1"]["baseMVA"]
                        cont_string = join(["tie line ", branch_id, " with ", cont_mw, " MW"])
                    elseif !isnothing(cont_name["dcbranch_id"])
                        dcbranch_id = cont_name["dcbranch_id"]
                        cont_mw = result_dc["solution"]["nw"]["$hour"]["branchdc"]["$dcbranch_id"]["pf"] * input_data["nw"]["1"]["baseMVA"]
                        cont_string = join(["dc line ", dcbranch_id, " with ", cont_mw, " MW"])
                    end
                    print("hour: ", idx, " contingency: ", cont_string, " dc contr = ", in_dc[idx, jdx], "\n")
                end
            end
        end
        in_dc_p = maximum(in_dc, dims = 2)
        f = join(["\$f_{min} = \$", fmin, " Hz"])
        if length(fmin) == 1
            Plots.plot!(p_in_hvdc, 1:length(hour_ids), in_dc_p[:, 1], linestyle = :dash, linecolor = :black, marker = :diamond, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$E^{dc}~in~GWs\$", label = f, xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern")
        else
            Plots.plot!(p_in_hvdc, 1:length(hour_ids), in_dc_p[:, 1], xlabel = "\$hour~id\$", ylabel = "\$E^{dc}~in~GWs\$", label = f, xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern")
        end
    end
    # p_in_hvdc = Plots.plot(1:length(hour_ids), in_dc_p[:, 1], linestyle = :dash, linecolor = :black, marker = :diamond, markercolor = :black, xlabel = "\$hour~id\$", ylabel = "\$E^{dc}~in~GWs\$", legend = false))
    plot_filename = joinpath("results", scenario, year, hours, join(["worst_case_hvdc_contributions.pdf"]))
    Plots.savefig(p_in_hvdc, plot_filename)
end

# function create_calc_time(fmin, scenario, year, hours)
#     t_opt_dc = Dict{String, Any}(["$i"=>[] for i in 1:length(fmin)])
#     t_opt_no_dc = Dict{String, Any}(["$i"=>[] for i in 1:length(fmin)])

#     for idx in 1:length(fmin)
#         fmin_ = fmin[idx]
#         print("Loading results", "\n")
#         fn = joinpath("results", scenario, year, hours, join(["f",fmin_,"_with_dc.json"]))
#         result_dc = Dict{String, Any}()
#         open(fn) do f
#         dicttxt = read(f,String)  # file information to string
#             result_dc = JSON.parse(dicttxt)  # parse and transform data
#         end

#         fn = joinpath("results", scenario, year, hours, join(["f",fmin_,"_without_dc.json"]))
#         result_no_dc = Dict{String, Any}()
#         open(fn) do f
#         dicttxt = read(f,String)  # file information to string
#             result_no_dc = JSON.parse(dicttxt)  # parse and transform data
#         end
#         t_opt_dc["$idx"] = result_dc["solve_time"]
#         t_opt_no_dc["$idx"] = result_no_dc["solve_time"]
#     end

#     filename = joinpath("results", scenario, year, hours, join(["solver_time_dc.json"]))
#     json_string = JSON.json(t_opt_dc)
#     open(filename,"w") do f
#     write(f, json_string)
#     end

#     filename = joinpath("results", scenario, year, hours, join(["solver_time_no_dc.json"]))
#     json_string = JSON.json(t_opt_no_dc)
#     open(filename,"w") do f
#     write(f, json_string)
#     end
# end


# p2 = Plots.plot(fmin, t_no_dc', marker = :diamond, xlabel = "\$f_{min} in~Hz\$", ylabel = "\$Calculation time~in~s\$", label = "without HVDC contribution")
# Plots.plot!(p2, fmin, t_dc', marker = :diamond,  xlabel = "\$f_{min} in~Hz\$", ylabel = "\$Calculation time~in~s\$", label = "with HVDC contribution")
# plot_filename = joinpath("results", "calculation_time_comparison.pdf")
# Plots.savefig(p2, plot_filename)



# fn = joinpath("results", join(["f",49.0,"_with_dc.json"]))
# result_dc = Dict{String, Any}()
# open(fn) do f
#     dicttxt = read(f,String)  # file information to string
#     global result_dc = JSON.parse(dicttxt)  # parse and transform data
# end

# fn = joinpath("results", join(["f",49.0,"_without_dc.json"]))
# result_no_dc = Dict{String, Any}()
# open(fn) do f
#     dicttxt = read(f,String)  # file information to string
#     global result_no_dc = JSON.parse(dicttxt)  # parse and transform data
# end

# number_of_generators = length(result_dc["solution"]["nw"]["1"]["gen"])
# number_of_hours = 24
# alpha_dc = zeros(number_of_hours, number_of_generators)
# alpha_p_dc = zeros(number_of_hours, number_of_generators)
# color_dc = Array{String}(undef, number_of_hours, number_of_generators)

# alpha_no_dc = zeros(number_of_hours, number_of_generators)
# alpha_p_no_dc = zeros(number_of_hours, number_of_generators)
# color_no_dc = Array{String}(undef, number_of_hours, number_of_generators)
# interval =  60 
# solutions =  1440
# gen_keys = sort(parse.(Int, collect(keys(result_dc["solution"]["nw"]["1"]["gen"]))))
# # for idx in 1:length(gen_keys)
# #     g = gen_keys[idx]
# #     alpha_dc[:, idx] = [result_dc["solution"]["nw"]["$idx"]["gen"]["$g"]["alpha_g"] for idx in 1:interval:solutions]
# #     alpha_no_dc[:, idx] = [result_no_dc["solution"]["nw"]["$idx"]["gen"]["$g"]["alpha_g"] for idx in 1:interval:solutions]
# # end

# # for i in 1:size(alpha_dc, 1)
# #     for j in 1:size(alpha_dc, 2)
# #         if alpha_dc[i,j] == 1
# #             color_dc[i,j] = "black"
# #         else
# #             color_dc[i,j] = "white"
# #         end
# #         if alpha_no_dc[i,j] == 1
# #             color_no_dc[i,j] = "black"
# #         else
# #             color_no_dc[i,j] = "white"
# #         end
# #         alpha_p_dc[i,j] = j
# #         alpha_p_no_dc[i,j] = j
# #     end
# # end

# # p3 = Plots.scatter(alpha_p_dc, color = color_dc, legend=nothing, xlabel = "\$hour\$", ylabel = "\$g_{id}\$", title = "generator status with HVDC contribution")
# # plot_filename = joinpath("results", "generator_status_dc.pdf")
# # Plots.savefig(p3, plot_filename)

# # p4 = Plots.scatter(alpha_p_no_dc,color=color_no_dc,legend=nothing, xlabel = "\$hour\$", ylabel = "\$g_{id}\$", title = "generator status without HVDC contribution")
# # plot_filename = joinpath("results", "generator_status_no_dc.pdf")
# # Plots.savefig(p4, plot_filename)


# # # for n in sort(parse.(Int, keys(result_dc["solution"]["nw"])))
# # #     print(n, " dc zone 1: ", result_dc["solution"]["nw"]["$n"]["zones"]["1"]["dc_contr"], " dc zone 2: ", result_dc["solution"]["nw"]["$n"]["zones"]["2"]["dc_contr"], "\n")
# # # end


# ls_dc = sum([sum([load["pcurt"] for (l, load) in result_dc["solution"]["nw"]["$idx"]["load"]]) for idx in 1:interval:solutions])
# ls_no_dc = sum([sum([load["pcurt"] for (l, load) in result_no_dc["solution"]["nw"]["$idx"]["load"]]) for idx in 1:interval:solutions])
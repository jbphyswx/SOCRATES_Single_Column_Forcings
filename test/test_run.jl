## Run test on cases and plot relative to LES

FT = Float64

flight_number = 1
forcing_type = :obs_data

if forcing_type == :obs_data
    forcing_str = "Obs"
elseif forcing_type == :ERA5_data
    forcing_str = "ERA5"
else
    error("forcing_type must be :obs_data or :ERA5_data")
end

case_name = "SOCRATES_RF" * string(flight_number, pad = 2) * "_" * forcing_str * "_data"

reload_environment = true
if reload_environment
    using Pkg
    Pkg.activate(expanduser("~/Research_Schneider/CliMa/TurbulenceConvection.jl/integration_tests/"))
    Pkg.develop(path = expanduser("~/Research_Schneider/CliMa/SOCRATESSingleColumnForcings.jl"))
    Pkg.activate(expanduser("~/Research_Schneider/CliMa/TurbulenceConvection.jl/"))
    Pkg.develop(path = expanduser("~/Research_Schneider/CliMa/SOCRATESSingleColumnForcings.jl"))
    Pkg.activate(expanduser("~/Research_Schneider/CliMa/TurbulenceConvection.jl/integration_tests/"))
    Pkg.develop(path = expanduser("~/Research_Schneider/CliMa/SOCRATESSingleColumnForcings.jl"))
    # Pkg.develop(path=expanduser("~/Research_Schneider/CliMa/TurbulenceConvection.jl/"))

    # Pkg.activate(expanduser("~/Research_Schneider/CliMa/SOCRATESSingleColumnForcings.jl/"))
    using TurbulenceConvection
    tc = pkgdir(TurbulenceConvection)
    include(joinpath(tc, "driver", "main.jl"))
    include(joinpath(tc, "driver", "generate_namelist.jl"))

end

setup_environment = true # leave true
if setup_environment

    namelist = NameList.default_namelist(
        case_name,
        root = "./",  # controls where the namelist gets written if write=true
        # write = parsed_args["write_namelist"], # controls if the namelist gets written (I'm not sure if this is repetitive -- does it get written again after the namelist overwrite? I think this is the wrong namelist that kept getting written when I ran things earlier separate from the results...)
        write = false, # don't write the namelist because we only care to write the later one once we've done our overwrites since that's the one that actually runs, though idk where that one gets written tbh...
    )

    namelist["time_stepping"]["dt_min"] = 0.5
    namelist["time_stepping"]["dt_max"] = 2.0
    namelist["time_stepping"]["t_max"] = 3600.0 * 1.1
    namelist["stats_io"]["frequency"] = 600.0
    supersat_type = :Base
    namelist["microphysics"]["τ_sub_dep"] = 10000.0
    namelist["microphysics"]["τ_cond_evap"] = 1.0
    namelist["user_args"] = (; use_supersat = supersat_type, τ_use = :morrison_milbrandt_2015_style)
    # namelist["user_args"] = (;use_supersat=supersat_type, τ_use=:morrison_milbrandt_2015_style_exponential_part_only) 
    # namelist["user_args"] = (;use_supersat=supersat_type, τ_use=:standard) 


    namelist["user_aux"] = Dict()
    namelist["user_aux"]["initial_profile_updraft_area"] = 0.33 # i think either dict or named tuple is fine, we dispatch depending...


    # namelist["turbulence"]["EDMF_PrognosticTKE"]["updraft_mixing_frac"] = 0.6
    # ("turbulence", "EDMF_PrognosticTKE", "max_area", FT(.3)), # stability limiting...
    # namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area_bc"] = "Prognostic" # unstable w/o  setting other params
    namelist["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = FT(0.7) # stability (maybe we need to use the limiter instead tho to not get flat cloud tops?)
    namelist["turbulence"]["EDMF_PrognosticTKE"]["min_area"] = FT(1e-10) # stability (maybe we need to use the limiter instead tho to not get flat cloud tops?)

    # namelist["turbulence"]["EDMF_PrognosticTKE"]["min_area"] = FT(0.3) # testing
    namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] = 0.35

    namelist["turbulence"]["EDMF_PrognosticTKE"]["limit_min_area"] = true # testing
    namelist["turbulence"]["EDMF_PrognosticTKE"]["min_area_limiter_scale"] = 5 # testing strong detrainment...
    namelist["turbulence"]["EDMF_PrognosticTKE"]["min_area_limiter_power"] = 2000 # testing strong detrainment...


    # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = 20 # testing strong detrainment...
    # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 4 # 

    namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = 20 # testing weak detrainment...
    namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 30 # 

    namelist["microphysics"]["τ_acnv_rai"] = 250.0
    namelist["microphysics"]["τ_acnv_sno"] = 250.0
    namelist["microphysics"]["q_liq_threshold"] = 1e-3
    namelist["microphysics"]["q_ice_threshold"] = 1e-5


    namelist["thermodynamics"]["moisture_model"] = "nonequilibrium"
    namelist["thermodynamics"]["sgs"] = "mean"

    namelist["output"]["output_root"] = expanduser("~/Research_Schneider/CliMa/SOCRATESSingleColumnForcings.jl/test")
    namelist["meta"]["uuid"] = "" # no uuid
    uuid = string(namelist["meta"]["uuid"])
    simname = namelist["meta"]["simname"]
    @info("simname", simname)
    outpath = joinpath(namelist["output"]["output_root"], "Output.$simname.$uuid") # no uuid
    @info("outpath", outpath)
    simulation_outpath = joinpath(outpath, "stats/Stats." * case_name * ".nc")
    # simulation_outpath = joinpath(outpath, "stats/Stats."*case_name*"_alt.nc") # testing an old one before any recent changes


    # output_relpath   =  "Output."*case_name*".out/stats/Stats."*case_name*".nc" # could try Glob.jl, see https://discourse.julialang.org/t/rdir-search-recursive-for-files-with-a-given-name-pattern/75605/20
    # simulation_outpath = joinpath(namelist["output"]["output_root"], output_relpath), 

    # @info(namelist)

end


run_simulation = true
if run_simulation
    main1d(namelist)
end



SOCRATES_truth_path = "/home/jbenjami/Research_Schneider/CliMa/SOCRATESSingleColumnForcings.jl/Data/Atlas_LES_Profiles/Output_Data/"
truth_file = joinpath(
    SOCRATES_truth_path,
    "RF" * string(flight_number, pad = 2) * "_" * forcing_str * "_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc",
)



# Read output file from SCM run and plot


using Dierckx
function pyinterp(x, xp, fp; bc = "error")
    spl = Dierckx.Spline1D(xp, fp; k = 1, bc = bc)
    return spl(vec(x))
end

using Plots
using NCDatasets
# create directory if doesn't exist
if !isdir(joinpath(outpath, "Figures"))
    mkdir(joinpath(outpath, "Figures"))
end

truth_data = NCDatasets.Dataset(truth_file, "r")
simul_data = NCDatasets.Dataset(simulation_outpath, "r")
z_simul = simul_data.group["profiles"]["zc"][:] # should be same for both but the truth one doesn't match the file atlas gave...
z_truth = truth_data["z"][:] # should be same for both but the truth one doesn't match the file atlas gave...
z = z_truth


# ============================================================================================================================================================= #
# create figure
simul_ind = 2
# ============================================================================================================================================================= #

t_simul = simul_data.group["timeseries"]["t"][:]
t_simul_ind = t_simul[simul_ind]

t_truth = truth_data["time"][:]
t_truth = (t_truth .- t_truth[1]) .* (24 * 3600) # make it relative to start time
truth_ind = argmin(abs.(t_truth .- t_simul_ind))
t_truth_ind = t_truth[truth_ind]



# ============================================================================================================================================================= #
# test box
# import CLIMAParameters as CP
# # import Thermodynamics as TD
# TD = TC.TD
# TDP = TC.TD.Parameters

# FT = Float64
# toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
# aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
# param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
# thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)

# Tg_simul_start = simul_data.group["timeseries"]["Tsurface"][1]
# Tg_simul_start = simul_data.group["profiles"]["temperature_mean"][:,1][1]
# pg = simul_data.group["reference"]["p_f"][:][1] # no timeseries
# ρg = TD.air_density(thermo_params,Tg_simul_start,pg)
# q_sfc = TD.q_vap_saturation_generic(thermo_params, Tg_simul_start, ρg, TD.Liquid()) # surface specific humidity
# q_sfc = TD.q_vap_saturation_from_pressure(thermo_params,0. , pg, Tg_simul_start, TD.PhaseEquil) # surface specific humidity








# ============================================================================================================================================================= #




# figure out simulation time we're plotting and what time that is in the truth_data (start should be same...)

### temperature ###
data_simul = simul_data.group["profiles"]["temperature_mean"][:, simul_ind]
data_truth = truth_data["TABS"][:, truth_ind]

ymax = 3000.0
ymin = -100.0


p_t = plot(
    data_truth,
    z_truth,
    label = "Truth",
    legend = :bottomleft,
    dpi = 600,
    marker = :circle,
    markersize = 0.5,
    markerstrokewidth = 0.2,
    xlabel = "T (K)",
    ylabel = "Height (m)",
    title = "Temperature (K)",
    ylim = (ymin, ymax),
)
p_s = plot!(data_simul, z_simul, label = "Simulation", marker = :circle, markersize = 0.5, markerstrokewidth = 0.2) # add simulation data to same plot
# vertical line at 5.22 + 273.15 degrees
# surface temp dot 
plot!(
    simul_data.group["timeseries"]["Tsurface"][[simul_ind]],
    [0],
    label = "Surface Temp",
    seriestype = :scatter,
    markersize = 5,
    markercolor = :red,
)
plot!(truth_data["SST"][[truth_ind]], [0], label = "SST", seriestype = :scatter, markersize = 5, markercolor = :blue)

savefig(joinpath(outpath, "Figures", "temperature.png")) # save to file

p_t = plot(
    data_simul .- pyinterp(z_simul, z_truth, data_truth, bc = "extrapolate"),
    z_simul,
    label = "Truth",
    legend = :bottomright,
    ylim = (ymin, ymax),
    xlabel = "T diff (K)",
    ylabel = "Height (m)",
    title = "Temperature (K) diff",
)
savefig(joinpath(outpath, "Figures", "temperature_diff.png")) # save to file

### qt_mean ###
data_simul = simul_data.group["profiles"]["qt_mean"][:, simul_ind]
data_truth = truth_data["QT"][:, truth_ind] ./ 1000

p_t = plot(
    data_truth,
    z_truth,
    label = "Truth",
    legend = :bottomleft,
    marker = :circle,
    markersize = 1.5,
    dpi = 600,
    ylim = (ymin, ymax),
    xlabel = "qt (kg/kg)",
    ylabel = "Height (m)",
    title = "QT (kg/kg)",
)
p_s = plot!(data_simul, z_simul, label = "Simulation", marker = :circle, markersize = 1.5) # add simulation data to same plot
savefig(joinpath(outpath, "Figures", "qt.png")) # save to file

p_t = plot(
    data_simul .- pyinterp(z_simul, z_truth, data_truth, bc = "extrapolate"),
    z_simul,
    label = "diff",
    legend = :bottomright,
    ylim = (ymin, ymax),
    xlabel = "qt diff (K)",
    ylabel = "Height (m)",
    title = "qt_diff",
)
savefig(joinpath(outpath, "Figures", "qt_diff.png")) # save to file



### ql_mean ### (need to figure out indexing here in particular..., interpolate in t)
ql_simul_ind = simul_ind
ql_truth_ind = truth_ind
# ql_ind = 1
data_simul = simul_data.group["profiles"]["ql_mean"][:, ql_simul_ind]
data_truth = truth_data["QCL"][:, ql_truth_ind] ./ 1000

p_t = plot(
    data_truth,
    z_truth,
    label = "Truth",
    legend = :bottomright,
    ylim = (ymin, ymax),
    xlabel = "ql (kg/kg)",
    ylabel = "Height (m)",
    title = "QL (kg/kg)",
)
p_s = plot!(data_simul, z_simul, label = "Simulation") # add simulation data to same plot
savefig(joinpath(outpath, "Figures", "ql.png")) # save to file

p_t = plot(
    data_simul .- pyinterp(z_simul, z_truth, data_truth, bc = "extrapolate"),
    z_simul,
    label = "diff",
    legend = :bottomright,
    ylim = (ymin, ymax),
    xlabel = "ql diff (K)",
    ylabel = "Height (m)",
    title = "ql_diff",
)
savefig(joinpath(outpath, "Figures", "ql_diff.png")) # save to file


### qi_mean ### (need to figure out indexing here in particular..., interpolate in t)
qi_simul_ind = simul_ind
qi_truth_ind = truth_ind
# ql_ind = 1
data_simul = simul_data.group["profiles"]["qi_mean"][:, qi_simul_ind]
data_truth = truth_data["QCI"][:, qi_truth_ind] ./ 1000

p_t = plot(
    data_truth,
    z_truth,
    label = "Truth",
    legend = :bottomright,
    ylim = (ymin, ymax),
    xlabel = "qi (kg/kg)",
    ylabel = "Height (m)",
    title = "QI (kg/kg)",
)
p_s = plot!(data_simul, z_simul, label = "Simulation") # add simulation data to same plot
savefig(joinpath(outpath, "Figures", "qi.png")) # save to file

p_t = plot(
    data_simul .- pyinterp(z_simul, z_truth, data_truth; bc = "extrapolate"),
    z_simul,
    label = "diff",
    legend = :bottomright,
    ylim = (ymin, ymax),
    xlabel = "qi diff (K)",
    ylabel = "Height (m)",
    title = "qi_diff",
)
savefig(joinpath(outpath, "Figures", "qi_diff.png")) # save to file


### qc_mean ### (need to figure out indexing here in particular..., interpolate in t)
qc_simul_ind = simul_ind
qc_truth_ind = truth_ind
# ql_ind = 1
data_simul =
    simul_data.group["profiles"]["ql_mean"][:, qc_simul_ind] + simul_data.group["profiles"]["qi_mean"][:, qc_simul_ind]
data_truth = truth_data["QCL"][:, qc_truth_ind] ./ 1000 + truth_data["QCI"][:, qc_truth_ind] ./ 1000

p_t = plot(
    data_truth,
    z_truth,
    label = "Truth",
    legend = :bottomright,
    ylim = (ymin, ymax),
    xlabel = "qc (kg/kg)",
    ylabel = "Height (m)",
    title = "QC (kg/kg)",
)
p_s = plot!(data_simul, z_simul, label = "Simulation") # add simulation data to same plot
savefig(joinpath(outpath, "Figures", "qc.png")) # save to file

p_t = plot(
    data_simul .- pyinterp(z_simul, z_truth, data_truth; bc = "extrapolate"),
    z_simul,
    label = "diff",
    legend = :bottomright,
    ylim = (ymin, ymax),
    xlabel = "qc diff (K)",
    ylabel = "Height (m)",
    title = "qc_diff",
)
savefig(joinpath(outpath, "Figures", "qc_diff.png")) # save to file

# liquid potential Temp

data_simul = simul_data.group["profiles"]["thetal_mean"][:, simul_ind]
data_truth = truth_data["THETAL"][:, truth_ind]

p_t = plot(
    data_truth,
    z_truth,
    label = "Truth",
    legend = :bottomright,
    ylim = (ymin, ymax),
    xlabel = "theta_l (K)",
    ylabel = "Height (m)",
    title = "theta_l (K)",
)
p_s = plot!(data_simul, z_simul, label = "Simulation") # add simulation data to same plot
savefig(joinpath(outpath, "Figures", "theta_l.png")) # save to file

p_t = plot(
    data_simul .- pyinterp(z_simul, z_truth, data_truth; bc = "extrapolate"),
    z_simul,
    label = "diff",
    legend = :bottomright,
    ylim = (ymin, ymax),
    xlabel = "theta_l diff (K)",
    ylabel = "Height (m)",
    title = "thetal_diff",
)
savefig(joinpath(outpath, "Figures", "theta_l_diff.png")) # save to file


# ql and qi combined on separate axes
data_simul_ql = simul_data.group["profiles"]["ql_mean"][:, simul_ind]
data_truth_ql = truth_data["QCL"][:, truth_ind] ./ 1000

data_simul_qi = simul_data.group["profiles"]["qi_mean"][:, simul_ind]
data_truth_qi = truth_data["QCI"][:, truth_ind] ./ 1000

p_t = plot(
    data_truth_ql,
    z_truth,
    label = "Truth liq",
    legend = :bottomright,
    ylim = (ymin, ymax),
    color = :green,
    xlabel = "q (kg/kg)",
    ylabel = "Height (m)",
    title = "QL and QI (kg/kg)",
)
p_s = plot!(data_simul_ql, z_simul, label = "Simulation liq", color = :lime) # add simulation data to same plot

# add ice data on same plot with twin x axis
xaxis2 = Plots.twiny()
p_s = plot!(xaxis2, data_truth_qi, z_truth, label = "Truth ice", ylim = (ymin, ymax), legend = :bottom, color = :blue)
p_s = plot!(xaxis2, data_simul_qi, z_simul, label = "Simulation ice", color = :cyan) # add simulation data to same plot
savefig(joinpath(outpath, "Figures", "ql_qi.png")) # save to file


# ql, qi updraft/env
data_simul_ql_up = simul_data.group["profiles"]["updraft_ql"][:, simul_ind]
data_simul_qi_up = simul_data.group["profiles"]["updraft_qi"][:, simul_ind]
data_simul_ql_env = simul_data.group["profiles"]["env_ql"][:, simul_ind]
data_simul_qi_env = simul_data.group["profiles"]["env_qi"][:, simul_ind]

plot(
    data_simul_ql_up,
    z_simul,
    label = "Updraft liq",
    legend = :bottomright,
    ylim = (ymin, ymax),
    color = :green,
    xlabel = "q (kg/kg)",
    ylabel = "Height (m)",
    title = "QL and QI (kg/kg)",
)
plot!(data_simul_ql_env, z_simul, label = "Environment liq", color = :green, linestyle = :dash) # add simulation data to same plot
plot!(data_simul_qi_up, z_simul, label = "Updraft ice", color = :blue) # add simulation data to same plot
plot!(data_simul_qi_env, z_simul, label = "Environment ice", color = :blue, linestyle = :dash) # add simulation data to same plot
savefig(joinpath(outpath, "Figures", "ql_qi_updraft_env.png")) # save to file



# updraft and environmental temperature
data_simul_up = simul_data.group["profiles"]["updraft_temperature"][:, simul_ind]
data_simul_env = simul_data.group["profiles"]["env_temperature"][:, simul_ind]
data_simul = simul_data.group["profiles"]["temperature_mean"][:, simul_ind]

data_truth = truth_data["TABS"][:, truth_ind]

p_t = plot(
    data_simul_up,
    z_simul,
    label = "Updraft",
    legend = :bottomleft,
    ylim = (ymin, ymax),
    color = :red,
    xlabel = "T (K)",
    ylabel = "Height (m)",
    title = "Updraft and Environmental Temperature (K)",
)
p_s = plot!(data_simul_env, z_simul, label = "Environment", color = :blue) # add simulation data to same plot
p_s = plot!(data_simul, z_simul, label = "Simulation", color = :brown, xlim = (250, 280)) # add simulation data to same plot

plot!(
    simul_data.group["timeseries"]["Tsurface"][[simul_ind]],
    [0],
    label = "Surface Temp",
    seriestype = :scatter,
    markersize = 5,
    markercolor = :red,
)
plot!(truth_data["SST"][[truth_ind]], [0], label = "SST", seriestype = :scatter, markersize = 2, markercolor = :blue)

plot!(data_truth, z_truth, label = "Truth mean", color = :green) # add simulation data to same plot

savefig(joinpath(outpath, "Figures", "updraft_env_temp.png")) # save to file




# updraft and environmental temperature start and now
data_simul_up = simul_data.group["profiles"]["updraft_temperature"][:, simul_ind]
data_simul_env = simul_data.group["profiles"]["env_temperature"][:, simul_ind]
data_truth = truth_data["TABS"][:, truth_ind]

data_simul_up_start = simul_data.group["profiles"]["updraft_temperature"][:, 1]
data_simul_env_start = simul_data.group["profiles"]["env_temperature"][:, 1]
data_truth_start = truth_data["TABS"][:, 1]



p_t = plot(
    data_simul_up,
    z_simul,
    label = "Updraft",
    legend = :bottomleft,
    ylim = (ymin, ymax),
    color = :red,
    dpi = 600,
    xlabel = "T (K)",
    ylabel = "Height (m)",
    title = "Updraft and Environmental Temperature (K)",
)
p_s = plot!(data_simul_env, z_simul, label = "Environment", color = :blue) # add simulation data to same plot

plot!(
    simul_data.group["timeseries"]["Tsurface"][[simul_ind]],
    [0],
    label = "Surface Temp",
    seriestype = :scatter,
    markersize = 5,
    markercolor = :red,
)
plot!(truth_data["SST"][[truth_ind]], [0], label = "SST", seriestype = :scatter, markersize = 2, markercolor = :blue)
plot!(data_truth, z_truth, label = "Truth mean", color = :green) # add simulation data to same plot

plot!(data_simul_up_start, z_simul, label = "Updraft start", color = :pink, linestyle = :dash) # add simulation data to same plot
plot!(data_simul_env_start, z_simul, label = "Environment start", color = :cyan, linestyle = :dash) # add simulation data to same plot
plot!(data_truth_start, z_truth, label = "Truth start", color = :lime, linestyle = :dash, xlim = (250, 280)) # add simulation data to same plot

savefig(joinpath(outpath, "Figures", "updraft_env_temp_start_now.png")) # save to file


# qt start vs now
data_simul = simul_data.group["profiles"]["qt_mean"][:, simul_ind]
data_simul_plus_precip =
    simul_data.group["profiles"]["qt_mean"][:, simul_ind] +
    simul_data.group["profiles"]["qs_mean"][:, simul_ind] +
    simul_data.group["profiles"]["qr_mean"][:, simul_ind]

data_simul_up = simul_data.group["profiles"]["updraft_qt"][:, simul_ind]
data_simul_env = simul_data.group["profiles"]["env_qt"][:, simul_ind]
data_truth =
    truth_data["QT"][:, truth_ind] ./ 1000 +
    truth_data["QS"][:, truth_ind] ./ 1000 +
    truth_data["QR"][:, truth_ind] ./ 1000

data_simul_start = simul_data.group["profiles"]["qt_mean"][:, 1]
data_simul_up_start = simul_data.group["profiles"]["updraft_qt"][:, 1]
data_simul_env_start = simul_data.group["profiles"]["env_qt"][:, 1]
data_truth_start = truth_data["QT"][:, 1] ./ 1000


p_t = plot(
    data_truth,
    z_truth,
    label = "Truth",
    legend = :bottomleft,
    ylim = (ymin, ymax),
    dpi = 600,
    marker = :circle,
    markersize = 1.5,
    xlabel = "qt (kg/kg)",
    ylabel = "Height (m)",
    title = "QT (kg/kg)",
)
p_s = plot!(data_simul, z_simul, label = "Simulation", marker = :circle, markersize = 1.5) # add simulation data to same plot
plot!(data_simul_up, z_simul, label = "Updraft", color = :red) # add simulation data to same plot
plot!(data_simul_env, z_simul, label = "Environment", color = :blue) # add simulation data to same plot
plot!(data_simul_plus_precip, z_simul, label = "simul+Precip", color = :black) # add simulation data to same plot

plot!(data_truth_start, z_truth, label = "Truth start", color = :green, linestyle = :dash) # add simulation data to same plot
plot!(data_simul_start, z_simul, label = "Simulation start", color = :pink, linestyle = :dash) # add simulation data to same plot
plot!(data_simul_up_start, z_simul, label = "Updraft start", color = :red, linestyle = :dash) # add simulation data to same plot
plot!(data_simul_env_start, z_simul, label = "Environment start", color = :blue, linestyle = :dash) # add simulation data to same plot

savefig(joinpath(outpath, "Figures", "qt_start_now.png")) # save to file


# qt start vs now vs pressure
# data_simul = simul_data.group["profiles"]["qt_mean"][:,simul_ind] 
# data_simul_plus_precip = simul_data.group["profiles"]["qt_mean"][:,simul_ind]  + simul_data.group["profiles"]["qs_mean"][:,simul_ind]  + simul_data.group["profiles"]["qr_mean"][:,simul_ind]

# data_simul_up = simul_data.group["profiles"]["updraft_qt"][:,simul_ind]
# data_simul_env = simul_data.group["profiles"]["env_qt"][:,simul_ind]
data_truth = truth_data["QT"][:, truth_ind] ./ 1000

# data_simul_start = simul_data.group["profiles"]["qt_mean"][:,1]
# data_simul_up_start = simul_data.group["profiles"]["updraft_qt"][:,1]
# data_simul_env_start = simul_data.group["profiles"]["env_qt"][:,1]
data_truth_start = truth_data["QT"][:, 1] ./ 1000

p_truth_start = truth_data["PRES"][:, 1]

p_t = plot(
    data_truth_start,
    p_truth_start,
    label = "Truth",
    legend = :bottomleft,
    dpi = 600,
    marker = :circle,
    markersize = 1.5,
    yflip = true,
    xlabel = "qt (kg/kg)",
    ylabel = "Height (m)",
    title = "QT (kg/kg)",
)
p_s = plot!(data_truth, p_truth_start, label = "Truth now", marker = :circle, markersize = 1.5)


savefig(joinpath(outpath, "Figures", "qt_start_now_vs_p.png")) # save to file




# updraft area
data_simul = simul_data.group["profiles"]["updraft_area"][:, simul_ind]

plot(
    data_simul,
    z_simul,
    label = "Updraft Area",
    legend = :topright,
    ylim = (ymin, ymax),
    xlabel = "Area",
    ylabel = "Height (m)",
    title = "Updraft Area",
)
savefig(joinpath(outpath, "Figures", "updraft_area.png")) # save to file

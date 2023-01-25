param_set = SOCRATES_Single_Column_Forcings.Parameters.param_set

function process_case(flight_number::Int; obs_or_ERA5 = "Obs"::Union{String,Symbol}, new_z::Union{Nothing,AbstractArray}=nothing, initial_condition::Bool=false, param_set=param_set, surface=nothing)
    """
    Processes the flight data for that case
    If new_z is not specified, will use Rachel Atlas default grid
    If initial_condition is false, returns the data at time 0, if true, returns full spline interpolations over the time dimension

    Some things are always forced by ERA5, otherwise we pull from whatever obs_or_ERA5 specifies... (list which is which?)
    """

    # initial conditions
    data  = SOCRATES_Single_Column_Forcings.open_atlas_les_profile(flight_number);

    if isnothing(new_z)
        new_z = data[:grid_data]
    end

    data = data[(:obs_data,:ERA5_data)]

    if obs_or_ERA5 == "Obs"
        forcing = :obs_data
    elseif obs_or_ERA5 == "ERA5"
        forcing = :ERA5_data
    else
        if obs_or_ERA5 ∉ [:obs_data, :ERA5_data, "Obs","ERA5"]
            error("obs_or_ERA5 must be either \"Obs\"/[:obs_data] or \"ERA5\"/[:ERA5_data]")
        else
            forcing = obs_or_ERA5
        end
    end

    # @show(obs_or_ERA5)

    # specify our action dimensions
    z_dim_num    = get_dim_num( "lev", data[forcing]["T"])
    time_dim_num = get_dim_num("time", data[forcing]["T"])


    
    # We use map here on individual variables mostly as many of them are no longer constituent variables of data so we must construct fcns with them separately

    data = data[(:obs_data,:ERA5_data)]


    if ~isnothing(surface)
        if surface ∈ ["reference_state", "reference","ref"]  # we just want the surface reference state and we'll just return that
            Tg = data[forcing]["Tg"][:][1] # might have to drop lon,lat dims or sum
            pg = data[forcing]["Ps"][:][1]
            qg = calc_qg(Tg, pg)
            return TD.PhaseEquil_pTq(thermo_params, pg , Tg , qg )
        elseif surface ∈ ["surface_conditions", "conditions","cond"] 
            Tg = vec(data[forcing]["Tg"])[:] # might have to drop lon,lat dims or sum
            pg = vec(data[forcing]["Ps"])[:]
            qg = calc_qg(Tg, pg)
            return (;pg=t->pyinterp(t, data[forcing]["tsec"][:], pg), Tg=t->pyinterp(t, data[forcing]["tsec"][:], Tg), qg=t->pyinterp(t, data[forcing]["tsec"][:], qg) )
                   
        else
            error("if surface is not set to nothing (i.e you want a surface value, it must be either \"reference_state\" or \"surface_conditions\"")
        end
    end
    
    # p,T,q | need these from both ERA and obs to construct ts, tsg
    p  = map(x->x["lev"], data)
    p  = map(x->align_along_dimension(x, z_dim_num), p) # align on dimension 3 lev (hopefully order was already correct (testing, can't recall if theres anywhere i didn't want this done...)
    pg = map(x->x["Ps"],  data)
    T  = map(x->x["T"], data)
    Tg = map(x->x["Tg"], data)
    q  = map(x->x["q"], data)


    qg       = map(calc_qg, Tg, pg)
    q_full   = combine_air_and_ground_data(q[forcing], qg[forcing], z_dim_num) # full q_array/qt_nudge -- is used from the chosen forcing dataset
    qt_nudge = q_full 

    # Set up thermodynamic states for easier use (for both forcing and ERA -- note ERA subsidence for example depends on density which relies on T,p,q so need both even if forcing is :obs_data)
    ts      = map((p,T,q)->TD.PhaseEquil_pTq.(thermo_params, p , T , q ), p,T,q)
    tsg     = map((pg,Tg,qg)->TD.PhaseEquil_pTq.(thermo_params, pg , Tg , qg ), pg,Tg,qg)
    ts_full = map((ts,tsg)->combine_air_and_ground_data(ts, tsg, z_dim_num), ts,tsg)

    # old_z  => Precompute old z coordinate (precompute to save us some trouble later (get_data_new_z_t func can self-calculate it but it's redundant to keep calculating z)
    z_old = map((ts,tsg,data)->lev_to_z( ts,tsg; param_set=param_set, data=data) , ts,tsg, data)

    # ω (subsidence) # always forced by era5
    ρ  = TD.air_density.(thermo_params, ts[:ERA5_data])
    ρg = TD.air_density.(thermo_params,tsg[:ERA5_data])
    ρ  = combine_air_and_ground_data(ρ, ρg, z_dim_num)
    ω  = data[:ERA5_data]["omega"] 
    dpdt_g = data[:ERA5_data]["Ptend"]
    dpdt_g = add_dim(dpdt_g, z_dim_num) # should be lon lat lev time (hopefully order was already correct)
    ω  = combine_air_and_ground_data(ω, dpdt_g, z_dim_num)

    p_grid = add_dim(align_along_dimension(p[:ERA5_data], z_dim_num),time_dim_num) # align on dimension 3 lev, add time dimensino 4
    # p_full = combine_air_and_ground_data(p_grid, add_dim(pg[:ERA5_data],z_dim_num) ,z_dim_num) # align p along lev, add lev_dim to pg, stack (hope braodcasting works to fill it out...)

    p_full = combine_air_and_ground_data(p_grid, pg[:ERA5_data][:] ,z_dim_num) # align p along lev, add lev_dim to pg, stack (hope braodcasting works to fill it out...)

    c    = 100
    a    = -2 + exp(250/c)
    f_p  = @. 2(a+1) / (a+exp(p_full/c)) - 1
    grav = TCP.grav(param_set)  # TD.Parameters.grav(thermo_params)
    subsidence =  -(ω .- dpdt_g.*f_p)  ./ (ρ .* grav) 
    
    # u, v # always forced by ERA5
    u = data[:ERA5_data]["u"] 
    v = data[:ERA5_data]["v"]
    u = combine_air_and_ground_data(u,FT(0),z_dim_num)
    v = combine_air_and_ground_data(v,FT(0),z_dim_num)
    u_nudge,v_nudge = u,v

    #=
    not sure what 
        "forced by ERA5-derived geostrophic winds and nudged toward the ERA5 horizontal winds with a nudging
        timescale of 1 hr for the ERA5-based simulation and 20 min for the Obs-based simulation."
    means so we're just using u,v not ug,vg for now...
    What does it mean to have both the forcing and the nudging? are they not overlapping/competing?
    =#

    # # u_g, v_g # always forced by ERA5 
    # ug = data[:ERA5_data]["ug"] 
    # vg = data[:ERA5_data]["vg"]
    # ug = combine_air_and_ground_data(ug,FT(0),z_dim_num)
    # vg = combine_air_and_ground_data(vg,FT(0),z_dim_num)

    # H (nudge)
    θ_liq_ice  = TD.liquid_ice_pottemp.(thermo_params, ts_full[forcing])
    H_nudge    = θ_liq_ice


    # dTdt_hadv (i think these are supposed to be always ERA5)
    dTdt_hadv = data[:ERA5_data]["divT"] 
    dTdt_hadv = combine_air_and_ground_data(dTdt_hadv,FT(0),z_dim_num)

    # dqdt_hadv (i think these are supposed to be always ERA5)
    dqtdt_hadv = data[:ERA5_data]["divq"] 
    dqtdt_hadv = combine_air_and_ground_data(dqtdt_hadv,FT(0),z_dim_num)

    # expand everything to new z grid and make t operation (initial condition or time splines)
    dTdt_hadv  = get_data_new_z_t(dTdt_hadv , new_z, z_dim_num,time_dim_num; z_old = z_old[:ERA5_data], data=data[:ERA5_data], param_set=param_set,  initial_condition=initial_condition)
    H_nudge    = get_data_new_z_t(H_nudge   , new_z, z_dim_num,time_dim_num; z_old = z_old[forcing]   , data=data[forcing]   , param_set=param_set,  initial_condition=initial_condition)
    dqtdt_hadv = get_data_new_z_t(dqtdt_hadv, new_z, z_dim_num,time_dim_num; z_old = z_old[:ERA5_data], data=data[:ERA5_data], param_set=param_set,  initial_condition=initial_condition)
    qt_nudge   = get_data_new_z_t(qt_nudge  , new_z, z_dim_num,time_dim_num; z_old = z_old[forcing]   , data=data[forcing]   , param_set=param_set,  initial_condition=initial_condition)
    subsidence = get_data_new_z_t(subsidence, new_z, z_dim_num,time_dim_num; z_old = z_old[:ERA5_data], data=data[:ERA5_data], param_set=param_set,  initial_condition=initial_condition)
    u_nudge    = get_data_new_z_t(u_nudge   , new_z, z_dim_num,time_dim_num; z_old = z_old[:ERA5_data], data=data[:ERA5_data], param_set=param_set,  initial_condition=initial_condition)
    v_nudge    = get_data_new_z_t(v_nudge   , new_z, z_dim_num,time_dim_num; z_old = z_old[:ERA5_data], data=data[:ERA5_data], param_set=param_set,  initial_condition=initial_condition)


    return  (; dTdt_hadv, H_nudge, dqtdt_hadv, qt_nudge, subsidence, u_nudge, v_nudge)

end


function get_default_new_z(flight_number::Int;)
    data  = SOCRATES_Single_Column_Forcings.open_atlas_les_profile(flight_number);
    new_z = data[:grid_data]
    return new_z
end

function surface_ref_state(flight_number::Int; obs_or_ERA5 = "Obs"::Union{String,Symbol}, param_set=param_set)
    """
    In Cases.jl from TurbulenceConvection.jl you're asked to provide a surface reference state so I'm gonna provide one based off just whatever our base profile is
    """


end





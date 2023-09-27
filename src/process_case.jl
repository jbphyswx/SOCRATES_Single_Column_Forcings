
"""
    process_case(
        flight_number::Int;
        obs_or_ERA5 = "Obs"::Union{String,Symbol},
        new_z::Union{Nothing,AbstractArray}=nothing,
        initial_condition::Bool=false,
        thermo_params,
        surface=nothing
    )

Processes the flight data for that case
If new_z is not specified, will use Rachel Atlas default grid
If initial_condition is false, returns the data at time 0, if true, returns full spline interpolations over the time dimension

Some things are always forced by ERA5, otherwise we pull from whatever obs_or_ERA5 specifies... (list which is which?)
"""
function process_case(
        flight_number::Int;
        obs_or_ERA5 = "Obs"::Union{String,Symbol},
        new_z::Union{Nothing,AbstractArray}=nothing,
        initial_condition::Bool=false,
        thermo_params,
        surface=nothing
    )

    ##### ALL OUR combine_air_and_ground_data calls need to be retrofitted to actually account for where the ground is! ######

    ## Plan -- create a stored variable "combine_z_dim" or something like that that you can pass around that tells you where to insert the values from ground into the full array in combine_air_and_ground data -- and some argument that uses that variable to supplant other options (or make another func?)

    # initial conditions
    data  = open_atlas_les_input(flight_number);

    if isnothing(new_z)
        new_z = data[:grid_data]
    end


    if obs_or_ERA5 ∈ ["Obs", :obs_data]
        forcing = :obs_data
        data = data[(:obs_data,:ERA5_data)]
    elseif obs_or_ERA5 ∈ ["ERA5", :ERA5_data]
        forcing = :ERA5_data
        data = data[(:ERA5_data,)] # drop obs if we're doing era5 cause we don't need it (11 has no obs either i think)        
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

    ## For surface values, adjust to the time period of interest and return only the surface values, these are called in TC.jl
    if ~isnothing(surface) #(is always ERA5) -- maybe not, doesn't seem so in ATLAS LES outputs
        initial_ind = get_initial_ind(data[forcing], flight_number) # should the reference be the first timestep? or what is the reference meant to be...
        if surface ∈ ["reference_state", "reference","ref"]  # we just want the surface reference state and we'll just return that
            Tg = data[forcing]["Tg"][:][initial_ind] # might have to drop lon,lat dims or sum
            pg = data[forcing]["Ps"][:][initial_ind]
            # qg = calc_qg(Tg, pg; thermo_params)

            p = vec(data[forcing]["lev"])[:]
            q = vec(selectdim(data[forcing]["q"], time_dim_num, initial_ind))[:] # select our q value subset along the time dimension

            qg = calc_qg([pg],p,q)
            qg = collect(qg)[]  # Thermodynamics 0.10.2 returns a tuple rather than scalar, so this can collapse to scalar in either 0.10.1<= or 0.10.2>=
            return TD.PhaseEquil_pTq(thermo_params, pg , Tg , qg )
        elseif surface ∈ ["surface_conditions", "conditions","cond"]
            # Tg = vec(data[forcing]["Tg"])[:][initial_ind:end] # might have to drop lon,lat dims or sum
            pg = vec(data[forcing]["Ps"])[:][initial_ind:end]
            # qg = calc_qg(Tg, pg; thermo_params)

            p = vec(data[forcing]["lev"])[:]
            q = selectdim(data[forcing]["q"], time_dim_num, initial_ind:size(data[forcing]["q"], time_dim_num)) # select our q value subset along the time dimension
            q = vec.(collect(eachslice(q, dims=time_dim_num))) # turn our q from [lon, lat, lev, time] to a list of vectors along [lon,lat,lev] to match p
            qg = map((pg,q)->calc_qg([pg],p,q)[1], pg, q) # map the function to get out qg for each time step

            tg = data[forcing]["tsec"][initial_ind:end] # get the time array
            tg = tg .- tg[1] # i think we need this to get the initial time to be 0, so the interpolation works
            # in this interpolation, tsec has to be adjusted to our offsets no? or we clip e.g. pg to be pg[initial_ind:end], tsec would also need to be adjusted no?, subtract the value at initial_ind i guess...
            return (;pg=t->pyinterp([t], tg, pg)[1], Tg=t->pyinterp([t], tg, Tg)[1], qg=t->pyinterp([t], tg, qg)[1] ) # would use ref and broadcast but doesnt convert back to array
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
    # qg = map((Tg,pg)->calc_qg(Tg,pg;thermo_params), Tg, pg)
    qg = map((pg,p,q)->calc_qg(pg,p,q), pg, p, q)


    # Set up thermodynamic states for easier use (for both forcing and ERA -- note ERA subsidence for example depends on density which relies on T,p,q so need both even if forcing is :obs_data)
    # ts      = map((p,T,q)->TD.PhaseEquil_pTq.(thermo_params, p , T , q ), p,T,q)
    tsg     = map((pg,Tg,qg)->TD.PhaseEquil_pTq.(thermo_params, pg , Tg , qg ), pg,Tg,qg)

    #=
     I believe we can get from T_L to something better
    In principle at initiation they assume domination by liquid, so T_l,i ≈ T_l
    We are given in the forcing files T_L = T - L q_c/ c_p, 
    Then, θ_L,I ≈ θ_L = θ - θ\T L q_c/ c_p = (θ/T) ( T - L q_c/ c_p) = (θ/T) T_L
    thus θ_L = T_L (p_0/p)^k, where k = R_d/c_p
    In principle this is the same as just calculating the dry potential temperature subsitututing T_L for t (not sure if to include q in the calculation of θ or no, in priciple it's needed for R\c_p so I'll do it as a phase partition of all vapor right now)
    I think for the surface we're fine because they gave us real temperature...
    This I hope solves the issue where our qt is right but q_l and T are wrong
    =#
    θ  = map((T,p,q)->TD.dry_pottemp_given_pressure.(thermo_params,T,p,TD.PhasePartition.(q)), T,p,q)
    ts = map((p,θ,q)->TD.PhaseEquil_pθq.(thermo_params, p , θ, q ), p,θ,q)
    # tsg     = map((pg,Tg,qg)->TD.PhaseEquil_pθq.(thermo_params, pg , Tg , qg ), pg,Tg,qg)


    # get indices where the ground would get inserted in... (and convert to same shape/dims.)
    ground_indices = map( (ts) -> get_ground_insertion_indices(ts, tsg[forcing], z_dim_num; thermo_params, data=data), ts) # always use ERA5 surface...
    # @show(ground_indices, map(size,ground_indices))

    q_full   = combine_air_and_ground_data(q[forcing], qg[forcing], z_dim_num; insert_location=ground_indices[forcing]) # full q_array/qt_nudge -- is used from the chosen forcing dataset
    qt_nudge = q_full


    ts_full = map((ts,ground_indices)->combine_air_and_ground_data(ts, tsg[forcing], z_dim_num; insert_location=ground_indices), ts, ground_indices) # is it era5? their output les plots sure don't look it...

    # old_z  => Precompute old z coordinate (precompute to save us some trouble later (get_data_new_z_t func can self-calculate it but it's redundant to keep calculating z)
    # z_old = map((ts,tsg,data)->lev_to_z( ts,tsg; thermo_params, data=data) , ts,tsg, data) # should this be tsg[:ERA5_data] cause surface is always ERA5
    z_old = map((ts,data)->lev_to_z( ts,tsg[forcing]; thermo_params, data=data) , ts, data) # should this be tsg[:ERA5_data] cause surface is always ERA5 (is it?)
    # @show(z_old)


    # ω (subsidence) # always forced by era5
    ρ  = TD.air_density.(thermo_params, ts[forcing])
    ρg = TD.air_density.(thermo_params,tsg[forcing])
    ρ  = combine_air_and_ground_data(ρ, ρg, z_dim_num; insert_location=ground_indices[forcing])
    ω  = data[forcing]["omega"]
    dpdt_g = data[forcing]["Ptend"]
    dpdt_g = add_dim(dpdt_g, z_dim_num) # should be lon lat lev time (hopefully order was already correct)
    ω  = combine_air_and_ground_data(ω, dpdt_g, z_dim_num; insert_location=ground_indices[forcing])

    p_grid = map((p,)->add_dim(align_along_dimension(p, z_dim_num),time_dim_num), p) # align on dimension 3 lev, add time dimensino 4
    p_full = map((p_grid, pg, ground_indices)->combine_air_and_ground_data(p_grid, pg[:] ,z_dim_num; insert_location=ground_indices), p_grid, pg, ground_indices) # align p along lev, add lev_dim to pg, stack (hope braodcasting works to fill it out...)

    # this i think was wrong -- maybe get atlas to confirm what version of the sigmoid she used.
    # c    = 100
    # a    = -2 + exp(250/c)
    # f_p  = @. 2(a+1) / (a+exp(p_full/c)) - 1

    # squeeze = (x) -> dropdims(x, dims = (findall(size(x) .== 1)...,))

    # this was derived for a scaled sigmoid passing through (ps,1), (25000 Pa, 0)
    L = 2.2 # maximum value (shape parameter)
    a = -L/2
    p0 = 250. * 100
    # @show( pg[:ERA5_data][:], size( pg[:ERA5_data][:]) )
    p1 = add_dim(  pg[forcing][:], z_dim_num)
    # @show(p1, size(p1))
    (p1,y1) = (p1 , 1) # use our ground pressure (can't use pfull cause would need to pull)
    (p2,y2) = (p0                , 0)
    k = @. log((L/2 + 1)/(L/2 - 1)) / (p1-p0) # ps is an array so we have an array of ks
    # @show(squeeze(k),size(k))
    # @show(size(p_full))
    f_p = @. a + L  / (1 + exp(-k*(p_full[forcing]-p0)))


    grav = TDP.grav(thermo_params)
    # @show(size(f_p), size(ω), size(ρ))
    subsidence =  -(ω .- dpdt_g.*f_p)  ./ (ρ .* grav)

    # @show(minimum(subsidence), maximum(subsidence))
    # subsidence = max.(subsidence,0)  # stability test (unstable) (all ascent)
    subsidence = min.(subsidence,0)  # stability test (stable!) ( all subsidence )
    # m = .007 # almost stable .005 | unstable .007
    # subsidence[subsidence .> m] .= m # stability test -- can we tolerate any ascent?

    # @show(squeeze(ω)', squeeze(f_p)' )
    # @show(squeeze(subsidence)',squeeze(p_full)')

    # u, v # always forced by ERA5
    u = data[forcing]["u"]
    v = data[forcing]["v"]
    u = combine_air_and_ground_data(u,FT(0),z_dim_num; insert_location=ground_indices[forcing])
    v = combine_air_and_ground_data(v,FT(0),z_dim_num; insert_location=ground_indices[forcing])
    u_nudge,v_nudge = u,v

    #=
    not sure what
        "forced by ERA5-derived geostrophic winds and nudged toward the ERA5 horizontal winds with a nudging
        timescale of 1 hr for the ERA5-based simulation and 20 min for the Obs-based simulation."
    means so we're just using u,v not ug,vg for now...
    What does it mean to have both the forcing and the nudging? are they not overlapping/competing?
    =#

    # u_g, v_g # always forced by ERA5
    ug = data[forcing]["ug"]
    vg = data[forcing]["vg"]
    ug = combine_air_and_ground_data(ug,FT(0),z_dim_num; insert_location=ground_indices[forcing])
    vg = combine_air_and_ground_data(vg,FT(0),z_dim_num; insert_location=ground_indices[forcing])
    ug_nudge,vg_nudge = ug,vg

    # H (nudge)
    θ_liq_ice  = TD.liquid_ice_pottemp.(thermo_params, ts_full[forcing])
    H_nudge    = θ_liq_ice


    # dTdt_hadv (i think these are supposed to be always ERA5)
    dTdt_hadv = data[forcing]["divT"]
    dTdt_hadv = combine_air_and_ground_data(dTdt_hadv,FT(0),z_dim_num; insert_location=ground_indices[forcing])

    # dqdt_hadv (i think these are supposed to be always ERA5)
    dqtdt_hadv = data[forcing]["divq"]
    dqtdt_hadv = combine_air_and_ground_data(dqtdt_hadv,FT(0),z_dim_num; insert_location=ground_indices[forcing])
    # expand everything to new z grid and make t operation (initial condition or time splines)
    # dTdt_hadv  = get_data_new_z_t(dTdt_hadv , new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing], data=data[forcing], thermo_params,  initial_condition)
    dTdt_hadv  = get_data_new_z_t(dTdt_hadv  .* (forcing==forcing) , new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing], data=data[forcing], thermo_params,  initial_condition)  # testing not using forcing for obs to see if fixes anything
    H_nudge    = get_data_new_z_t(H_nudge   , new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing]   , data=data[forcing]   , thermo_params,  initial_condition)
    # dqtdt_hadv = get_data_new_z_t(dqtdt_hadv, new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing], data=data[forcing], thermo_params,  initial_condition)
    dqtdt_hadv = get_data_new_z_t(dqtdt_hadv .* (forcing==forcing) , new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing], data=data[forcing], thermo_params,  initial_condition) # testing not using forcing for obs to see if fixes anything
    qt_nudge   = get_data_new_z_t(qt_nudge  , new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing]   , data=data[forcing]   , thermo_params,  initial_condition)
    subsidence = get_data_new_z_t(subsidence, new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing], data=data[forcing], thermo_params,  initial_condition)
    u_nudge    = get_data_new_z_t(u_nudge   , new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing], data=data[forcing], thermo_params,  initial_condition)
    v_nudge    = get_data_new_z_t(v_nudge   , new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing], data=data[forcing], thermo_params,  initial_condition)

    ug_nudge   = get_data_new_z_t(ug_nudge  , new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing], data=data[forcing], thermo_params,  initial_condition)
    vg_nudge   = get_data_new_z_t(vg_nudge  , new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing], data=data[forcing], thermo_params,  initial_condition)

    return  (; dTdt_hadv, H_nudge, dqtdt_hadv, qt_nudge, subsidence, u_nudge, v_nudge, ug_nudge, vg_nudge)

end


function get_default_new_z(flight_number::Int;)
    data  = open_atlas_les_input(flight_number);
    new_z = data[:grid_data]
    return new_z
end

"""
    surface_ref_state(
        flight_number::Int;
        obs_or_ERA5 = "Obs"::Union{String,Symbol},
        thermo_params
    )

In Cases.jl from TurbulenceConvection.jl you're asked to provide a surface reference state so I'm gonna provide one based off just whatever our base profile is
"""
function surface_ref_state(flight_number::Int; obs_or_ERA5 = "Obs"::Union{String,Symbol}, thermo_params)
    return nothing
end


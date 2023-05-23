#= ...
created jan 4 2023
@jbphyswx
... =#

##
# - Maybe swith to using DimensionalData or something so we can have derived quantities (i.e. ts <thermodynamic state>) as labeled arrays and not have to keep passing their parents around...
# -- seems to be a little cumbersome and no easy way to do netCDF reads and conversion to array/dataset like types
##

#= Simple linear interpolation function, wrapping Dierckx (copied from TC.jl)  =#
function pyinterp(x, xp, fp)
    spl = Dierckx.Spline1D(xp, fp; k = 1)
    return spl(vec(x))
end

"""
    swapaxes(a, dim1, dim2)

Swaps the two dimensions of an array.
"""
function swapaxes(a, dim1, dim2)
    perm = collect(1:ndims(a))
    perm[dim1] = dim2
    perm[dim2] = dim1
    return permutedims(a, perm)
end

"""
    align_along_dimension(v,dim)

Assumes v has one non-singleton dimension which will be aligned along dim.
if dim does not exist, v is expanded till dim exists
"""
function align_along_dimension(v,dim)

    sz = size(v)
    nd = ndims(v)
    if nd < dim # add dims if they don't exist
        v = reshape(v, sz..., fill(1, dim-nd)...)
    end

    # get nonsingleton dimension
    nonsingleton_ind = findall(sz .> 1)
    l_ns = length(nonsingleton_ind)
    if l_ns > 1
        error("v has more than one non-singleton dimension")
    elseif l_ns < 1
        error("v has no non-singleton dimension")
    end

    v = swapaxes(v, nonsingleton_ind[1], dim) # move the non-singleton dimension to the dim position
    return v
end

"""
    add_dim(vardata, dimnum)

Allows you to add a dimension to a variable at the specified location
"""
function add_dim(vardata, dimnum)
    # both data need to be labeled
    sz_vardata = collect(size(vardata)) # array
    insert!(sz_vardata,dimnum,1) # insert sz 1 at this location
    return reshape(vardata, sz_vardata...) # reshape
end

"""
    data_to_tsg(data; thermo_params)

Does this cause we don't have surface q values as of rn
"""
function data_to_tsg(data; thermo_params)

    Tg = data["Tg"]
    pg = data["Ps"]

    # not sure what to do at surface, so assuming saturation at surface
    pvg           = TD.saturation_vapor_pressure.(thermo_params, Tg, TD.Liquid())
    molmass_ratio = TD.molmass_ratio(thermo_params)

    qg = (1 / molmass_ratio) .* pvg ./ (pg .- pvg) #Total water mixing ratio at surface , assuming saturation

    tsg = TD.PhaseEquil_pTq.(thermo_params, pg, Tg, qg)
    return tsg
end

"""
    data_to_ts
"""
function data_to_ts(data; do_combine_air_and_ground_data=false, thermo_params) # add type data is ncdataset

    T = data["T"]
    q = data["q"]
    p = data["lev"]

    # put p/lev on correct dimension...
    p = align_along_dimension(p, get_dim_num("lev", data["T"]))

    ts = TD.PhaseEquil_pTq.(thermo_params, p, T, q)

    if do_combine_air_and_ground_data
        concat_dim = get_dim_num("lev", T) # phaseequil reduces us down to lev... can't seem to just apply ufunc...
        tsg = data_to_tsg(data; thermo_params)
        ts = combine_air_and_ground_data(ts,tsg, concat_dim)
    end

    return ts
end


function insert_sorted(vect, val; by=monotonic_checker, monotonic_rev=monotonic_rev)
    # @show(vect,val)
    index = searchsortedfirst(vect, val; by=monotonic_checker, rev=monotonic_rev) #find index at which to insert x
    return insert!(vect, index, val) #insert x at index
end

# function insert_sorted(vect, val;by=monotonic_checker,monotonic_rev=monotonic_rev)
#     index = searchsortedfirst(vect, val;  by=monotonic_checker, rev=monotonic_rev) #find index at which to insert x
#     return insert!(vect, index, val) #insert x at index
# end


"""
Because we might not have monotonic data, we need to be able to both insert the ground data into our lev array where applicable and calculate our dz accordingly...
I guess in principle you could throw out p < ps but then you wouldn't be able to return an even array out so either way this function would need to exist for padding and such...

tsz should be a one dimensional array consising of [ts..., tg]

# to do -- add capability to use precomputed indices (insert_location)
"""
function lev_to_z_column(tsz; thermo_params, data=data)

    ts  = tsz[1:end-1] # need to make it a vector I guess... (not sure if this screws wit output shape)
    tsg = tsz[end]
    R_d           = TDP.R_d(thermo_params)
    grav          = TDP.grav(thermo_params)

    # dimnames    = NC.dimnames(data["T"]) # use this as default cause calculating ts doesn't maintain dim labellings
    # lev_dim_num = findfirst(x->x=="lev",dimnames)
    L           = length(ts)

    # @show(size(ts),size(tsg))
    index = searchsortedfirst(ts, tsg;  by=x->TD.air_pressure(thermo_params,x), rev=false) # find where the ground value would be inserted...
    insert!(ts,index,tsg) # only seems to be an in place option...
    tsz       = ts # replace with the reordered version

    Tvz       = TD.virtual_temperature.(thermo_params, tsz) # virtual temp, we havent returned these for now...
    pz        = TD.air_pressure.(thermo_params, tsz)
    Lz        = L+1 # cause we extended it using the ground...
    Tv_bar    = Statistics.mean((Tvz[1:Lz-1], Tvz[2:Lz]))
    p_frac    = pz[2:Lz] ./ pz[1:Lz-1]
    dz        = @. (R_d * Tv_bar / grav) * log(p_frac)

    # sum up from the bottom then subtract the height of the ground
    z         = reverse(cumsum(reverse(dz))) #cumsum(dz) # grid is already defined from  (from Grid.jl)
    # s_sz      = collect(size(z)) # should now be same dims as T
    # s_sz[ldn] = 1 # we just want to add a single slice of zeros for the ground -- these are all ocean cases so this should be fine for now...
    z         = [z...,0] # is this the right order? seems so based on the cat below but idk... if so might have to flip index to be L - index or something like that? (seems to be so...)
    # after the padding, z for the ground should be at the right index, and we can just subtract it out from the array...
    z         = z .- z[index]
    return z
end

# convert pressure to altitude...
"""
    lev_to_z

TODO: document
"""
function lev_to_z( p::FT, T::FT, q::FT, pg::FT, Tg::FT, qg::FT; thermo_params, data) where {FT <: Real}
    ts            = TD.PhaseEquil_pTq.(thermo_params, p, T, q)
    tsg           = TD.PhaseEquil_pTq.(thermo_params, pg, Tg, qtg)
    return lev_to_z(ts, tsg; data=data, thermo_params)
end


"""
    lev_to_z

ts is thermodynamic state
tsg is thermodynamic state for ground

if assume monotonic, everything should already be in the right order and we can use the vectorized version, otherwise we will use lev_to_z column applied column by column w/ mapslices
"""
function lev_to_z(ts, tsg; thermo_params, data, assume_monotonic = false )

    dimnames    = NC.dimnames(data["T"]) # use this as default cause calculating ts doesn't maintain dim labellings
    lev_dim_num = findfirst(x->x=="lev",dimnames)
    ldn         = lev_dim_num
    L           = size(ts,lev_dim_num)

    if !assume_monotonic
        tsz = combine_air_and_ground_data(ts,tsg,ldn;data,reshape_ground=true, insert_location=:end) # here we just want to assume monotonic so we can pass those to this fcn
        z   = mapslices(x->lev_to_z_column(x;thermo_params,data), tsz; dims=ldn) # need to make a stack cause that's all mapslices can take...
    else

        R_d           = TDP.R_d(thermo_params) # TD.Parameters.R_d(thermo_params)
        grav          = TDP.grav(thermo_params)  # TD.Parameters.grav(thermo_params)


        # tsz       = cat(ts,tsg;dims=ldn)
        tsz       = combine_air_and_ground_data(ts,tsg,ldn;data=data,reshape_ground=true, insert_location=x->TD.air_pressure(thermo_params,x)) # this doesnt work cause sometimes ps is more than the lowest value in lev...

        #actually we need to split dz into pos or neg depending on whether or not it's above ground... maybe best just to have fcn that goes along each column...

        Tvz       = TD.virtual_temperature.(thermo_params, tsz) # virtual temp, we havent returned these for now...
        pz        = TD.air_pressure.(thermo_params, tsz)
        Lz        = L+1 # cause we extended it using the ground...
        Tv_bar    = Statistics.mean((selectdim(Tvz, ldn, 1:Lz-1), selectdim(Tvz, lev_dim_num, 2:Lz)))
        p_frac    = selectdim( pz, lev_dim_num, 2:Lz) ./ selectdim( pz, ldn, 1:Lz-1)
        dz        = @. (R_d * Tv_bar / grav) * log(p_frac)
        # @show(dz)

        z         = reverse(cumsum(reverse(dz;dims=ldn);dims=ldn);dims=ldn) #cumsum(dz) # grid is already defined from  (from Grid.jl)
        # @show(z)
        s_sz      = collect(size(z)) # should now be same dims as T
        s_sz[ldn] = 1 # we just want to add a single slice of zeros for the ground -- these are all ocean cases so this should be fine for now...
        _z0       = zeros(Float64,s_sz...)
        z         = cat(z, _z0; dims=ldn)
    end
    return z
end

function z_from_data(data; thermo_params)
    ts  = data_to_ts(data; thermo_params, do_combine_air_and_ground_data=false)
    tsg = data_to_tsg(data; thermo_params)
    return lev_to_z(ts, tsg; thermo_params, data)
end


"""
    get_ground_insertion_indices

Get the indices where the ground tsg would fit into the array ts...
"""
function get_ground_insertion_indices(ts,tsg, concat_dim; thermo_params, data=data)
    function mapslice_func(vect; thermo_params=thermo_params, by=x->TD.air_pressure(thermo_params,x))
        vardata  = vect[1:end-1]
        vardatag = vect[end]
        index = searchsortedfirst(vardata, vardatag; by=by, rev=false)
        return index
    end
    # reshape ( TODO!! # use add_dim )
    sz_tsg = collect(size(tsg)) # array
    insert!(sz_tsg,concat_dim,1) # insert sz 1 at this location
    tsg =  reshape(tsg, sz_tsg...) # reshape (just adds the singleton dimension in)
    # concat and calculate
    ts = cat(ts,tsg; dims=concat_dim)
    return mapslices(mapslice_func, ts; dims=[concat_dim])
end




"""
    combine_air_and_ground_data

var  : the data variable or string variable name for data in the air
varg : the data variable or string variable name for data on the ground
data : is a container from which the data can be accessed as a string in data[var(g)] form

reshape_ground : expand the ground data to the same size as bottom slice of the air data (so you can pass in for example a single scalar or similar reduced dimension varg)
assume_monotonic : assume that the ground value is actually below everything in the array... in reality sometimes we have for example a ground pressure above that of the minimum in this array... (should be faster than using insert_location )
# might not need this anymore cause if insert_location is integer that jut works

# atlas stated "We use hourly pressure level data interpolated onto a horizontal grid of 0.25° × 0.25° and 37 pressure levels from its native 137 hybrid sigma/pressure levels and 30 km horizontal grid." so i guess sometimes this causes slight problems
"""
function combine_air_and_ground_data(var,varg, concat_dim; data=nothing, reshape_ground=true, insert_location=:end)
    # if data is a string, read the data out from data (creates data `vardata` from `var` whether var is string or already is data)
    vardata    = isa(var ,String) ? data[var ] : var
    vardatag   = isa(varg,String) ? data[varg] : varg
    concat_dim = isa(concat_dim,String) ? get_dim_num(concat_dim, vardata) : concat_dim # string to numbers

    # handle the ground data alignment
    if reshape_ground
        if isa(vardatag, Number) # allow for you to pass a single number and still concat
            sz_vardatag = collect(size(vardata)) # array
            sz_vardatag[concat_dim] = 1
            vardatag = fill(vardatag, sz_vardatag... ) # create array full with just this one value
        elseif isa(vardatag, NC.CFVariable) # vardatag should not have lev as a dimension so we need to add it in (and I guess check the others are in the same order)
            dimnamesg = NC.dimnames(vardatag)
            dimnames  = NC.dimnames(vardata)
            # in same order as vardata just in case
            vardatag = reshape(vardatag, (size(vardatag)...,1)) # add trailing singleton for lev
            dimnamesg = [dimnamesg..., "lev"]
            vardatag = permutedims(vardatag, [findfirst(x->x==dim,dimnames) for dim in dimnamesg] ) # permute into the right order
        elseif isa(vardatag, AbstractArray) # an unlabeled array, i think we can't guarantee then that the lev axis exists at all or what existing dimensions are so this just assumes everything is correct except the lev dimension being there
            # assume we need to add a new dimension at the same location as in the full array and order otherwise is preserved (quick check looks ok)
            sz_vardatag = collect(size(vardatag)) # array
            if ndims(vardatag) == (ndims(vardata)-1)
                ## TODO USE add_dim
                insert!(sz_vardatag,concat_dim,1) # insert sz 1 at this location
                vardatag =  reshape(vardatag, sz_vardatag...) # reshape (just adds the singleton dimension in)
            elseif ndims(vardatag) != ndims(vardata)
                error("size mismath, var and varg should have the same number of dimensions or one less (for a missing lev dimension)")
            end
        end
    end

    # broadcast out to match sizes except for concat dim
    sz_vardata  = collect(size(vardata)) # array
    sz_vardatag = collect(size(vardatag)) # array
    num_repeat  = sz_vardatag .÷ sz_vardata # integer division
    num_repeatg = sz_vardata .÷ sz_vardata # integer division
    num_repeat[ concat_dim] = 1 # don't repeat along concat dim
    num_repeatg[concat_dim] = 1 # don't repeat along concat dim
    vardatag = repeat(vardatag, num_repeatg...)
    vardata  = repeat(vardata , num_repeat...)

    if insert_location == :end # here we assume the ground values are below the values we add automatically. (we used to use identity fcn but i think we don't need that)
        vardata = cat(vardata,vardatag;dims=concat_dim) # we concatenate it at the end...
    elseif isa(insert_location,Function) # use function to determine where to insert our ground values, uses search sorted assuming the original array is already sorted (speedup from binary search)
        mapslice_func = function(vect; by=insert_location) # write this way cause can't define func inside conditional unless anonymous?, see https://github.com/JuliaLang/julia/issues/15602 , https://stackoverflow.com/a/65660721
            # @show(vect)
            vardata  = vect[1:end-1]
            vardatag = vect[end]
            # @show(vardata, vardatag)
            return insert_sorted(vardata,vardatag; by=insert_location)
        end
        vardata = cat(vardata,vardatag;dims=concat_dim)
        vardata = mapslices(mapslice_func, vardata; dims=[concat_dim])

    elseif isa(insert_location, AbstractArray ) # use provided indices to determine where to input values...
        if isa(vardatag, Number) # single value, just splice in our value
            vardata = cat(selectdim(vardata, concat_dim, 1:insert_location-1),vardatag, selectdim(vardata, concat_dim, insert_location:size(vardata,concat_dim)); dims=concat_dim) # insert the slice there...
        elseif isa(vardatag, AbstractArray) # an unlabeled array, i think we can't guarantee then that the lev axis exists at all or what existing dimensions are so this just assumes everything is correct except the lev dimension being there
            # assume we need to add a new dimension at the same location as in the full array and order otherwise is preserved (quick check looks ok)
            mapslice_func = function(vect) # write this way cause can't define func inside conditional unless anonymous?, see https://github.com/JuliaLang/julia/issues/15602 , https://stackoverflow.com/a/65660721
                # @show(vect)
                vardata  = vect[1:end-2]
                vardatag = vect[end-1]
                insert_location = Int(vect[end]) # undo if was coerced to FT
                # insert (a little more complicated now cause we aren't exactly getting the same size output... dont think that actually matters for mapslices...)
                # @show(vardata, vardatag, insert_location)
                return insert!(vardata,insert_location, vardatag)
            end
            vardata = cat(vardata,vardatag,insert_location; dims=concat_dim)
            vardata = mapslices(mapslice_func, vardata; dims=[concat_dim])
        end
    else
        error("unsupported input type for variable insert_location") # catch what would otherwise silently fail and return vardata
    end

    return vardata
end




"""
    get_dim_num

get the dimension number of dim from nc_data
"""
function get_dim_num(dim,nc_data=nothing)

    # convert interp_dim to a number stored in interp_dim_num
    if isa(dim, Number) # already a string, just returns...
        dim_num = dim
    elseif isa(dim,String)
        if isa(nc_data, NC.CFVariable) # allow finding the number from an ncdataset and a string
            dimnames = NC.dimnames(nc_data)
            dim_num = findfirst(x->x==dim,dimnames)
        elseif isa(nc_data,NC.NCDataset)
            @warn("dimension number for a dataset is not well defined, use a speciifc NC.CFVariable instead")
        elseif isa(nc_data, AbstractArray)
            error("cannot find dimension $(dim) in unlabeled data nc_data, pass in a labeled NCDataset instead...")
        else
            error("unsupported input type for nc_data $(typeof(nc_data))")
        end
    end
    return dim_num
end



"""
    interp_along_dim

interpolation data

# var : the data to be interpolated
# interp_dim: the name or dimension number along which we will do interpolation

# interp_dim_in: the coordinate on which the data is currently aligned
# interp_dim_out: the coordinate on which we would like to interpolate the data to

the function will create splines along interp_dim_in...
- if we set interp_dim_out, we actually evaluate the function along interp_dim at locations interp_dim_out,otherwise, we just return our unevaluated spline functions

- data_func is a func to be applied to the raw data before it is processed, though perhaps it it most useful if interp_dim_out is unset and we are returning functions...
- vectorize_in means your input is an array and you need to loop over it too (as opposed to just being a fixed template vector)

To Do : decide types
"""
function interp_along_dim(var, interp_dim, interp_dim_in; interp_dim_out=nothing,data=nothing, data_func = nothing, interp_dim_in_is_full_array=true, reshape_ground=true, verbose=false)

    # if data is a string, read the data out from data (creates data `vardata` from `var` whether var is string or already is data)
    vardata    = isa(var ,String) ? data[var ] : var


    # get the interpolation dimension and combine air and ground data
    interp_dim_num = get_dim_num(interp_dim, vardata)
    # vardata        = combine_air_and_ground_data(vardata ,vardatag, interp_dim_num; data=nothing, reshape_ground=true) # combine air and ground data... (also resolves strings to data)

    if !isnothing(data_func) # apply data_func if we need to
        vardata  = data_func(vardata)
    end

    # mapslices to apply along timedim, see https://docs.julialang.org/en/v1/base/arrays/#Base.mapslices
    if !interp_dim_in_is_full_array
        if isnothing(interp_dim_out)
            return mapslices(d -> dd -> pyinterp(dd, interp_dim_in, d)      , vardata, dims=[interp_dim_num,]) # will return a lambda fcn that can be evaluated along that dimensoin
        else
            return mapslices(d -> pyinterp(interp_dim_out, interp_dim_in, d), vardata, dims=[interp_dim_num,]) # lambda fcn will evaluate
        end
    else # vectorize over input dim values as well as data (no support for vectorize over output dim yet)
        # stack on new catd dimension, then split apart inside the fcn call
        catd = ndims(vardata)+1
        _input = cat(interp_dim_in, vardata; dims=catd ) # although maybe the input is just a vector in which case this won't work..... we could just pass it in
        if isnothing(interp_dim_out)
            return dropdims( mapslices(d -> dd -> pyinterp(dd , d[:,1], d[:,2])     , _input, dims=[interp_dim_num,catd]); dims=catd) # wll return a lambda fcn that can be evaluated along that dimensoin
        else
            return dropdims( mapslices(d -> pyinterp(interp_dim_out, d[:,1], d[:,2]), _input, dims=[interp_dim_num,catd]); dims=catd) # lambda fcn will evaluate
        end
    end
end


##########################
#### -- note is this redundant? i think techincally we're create and instantateously evaluate the spline over and over -- maybe we want to just store the spline and then evaluate it at different points?
##########################



# function var_to_new_z(;var=var, varg=nothing, old_z=z, new_z=default_new_z, interp_dim=ldn, data_func=x->reverse(x;dims=ldn), kwargs...)
#     """
#     data interpolation to new z coordinate
#     """

#     interp_dim = isa(interp_dim, String) ? get_dim_num(interp_dim, data) : interp_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension

#     # make sure our data is numeric and combine air and ground if necessary...
#     vardata    = isa(var, String) ? data[var ] : var
#     if !isnothing(varg)
#         vardatag   = isa(varg,String) ? data[varg] : varg
#         vardata    = combine_air_and_ground_data(vardata ,vardatag, ldn; reshape_ground=true) # combine air and ground data... (also resolves strings to data)
#     end

#     return interp_along_dim(vardata, interp_dim, reverse(z; dims=interp_dim); interp_dim_out=new_z, data_func=data_func, kwargs...)
# end

"""
    var_to_new_coord

if coord_new is nothing, then will return functions...
"""
function var_to_new_coord(var,coord_in, interp_dim; coord_new=nothing, data=nothing, data_func=nothing, kwargs...)
    vardata    = isa(var ,String) ? data[var ] : var
    if ~isnothing(data)
        interp_dim = isa(interp_dim, String) ? get_dim_num(interp_dim, data) : interp_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
    else
        interp_dim = isa(interp_dim, String) ? get_dim_num(interp_dim, vardata) : interp_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
    end

    # evaluate interp_dim_in_is_full_array based on the size of the input... interp_dim_in is full array false is much faster cause dont have to double loop in vectorization...
    return interp_along_dim( vardata, interp_dim, coord_in; interp_dim_out=coord_new, data=data, data_func=data_func, interp_dim_in_is_full_array=(size(coord_in)==size(vardata)), kwargs...)
end


"""
    get_data_new_z_t

Take data from our base setup, interpolate it to new z, then create time splines based on the t we have....
to vectorize properly over z_new, it should be the same shape as vardata+vardata_g
"""
function get_data_new_z_t(var, z_new, z_dim, time_dim, flight_number; thermo_params, varg=nothing, z_old=nothing, t_old=nothing, data=nothing, initial_condition = false, assume_monotonic=false)

    # get the data and dimensions we're working on,
    vardata    = isa(var ,String) ? data[var ] : var
    if ~isnothing(varg)
        vardatag   = isa(varg,String) ? data[varg] : varg
    end

    # combine air and ground data
    if ~isnothing(data)
        z_dim_num = isa(z_dim, String) ? get_dim_num(z_dim, vardata) : z_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
        time_dim_num = isa(time_dim, String) ? get_dim_num(time_dim, vardata) : time_dim
    else
        z_dim_num = isa(z_dim, String) ? get_dim_num(z_dim, vardata) : z_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
        time_dim_num = isa(time_dim, String) ? get_dim_num(time_dim, vardata) : time_dim
    end


    if isnothing(z_old)
        z_old = z_from_data(data; thermo_params); # uses ground value to create the old z, pads bottom w/ 0
    end
    if isnothing(t_old)
        t_old = data["tsec"][:] # check this unit was right in the files (may need to make sure it's subtracting out the first timestep so starts at 0) -- do we need to align this on a dimension?
    end

    print(data["bdate"][:])
    t_base = Dates.DateTime(string(data["bdate"][:]), Dates.DateFormat("yymmdd")) + Dates.Year(2000) # the base Date (using bdate not nbdate cause nbdate seems to have  bug in flight 9 (extra 0 in month spot))
    t      = t_base .+ Dates.Second.(t_old) # the actual dates
    summary_file = joinpath(dirname(@__DIR__), "Data", "SOCRATES_summary.nc")
    SOCRATES_summary = NC.Dataset(summary_file,"r")
    flight_ind = findfirst(SOCRATES_summary["flight_number"][:] .== flight_number)
    initial_time = SOCRATES_summary["reference_time"][flight_ind] - Dates.Hour(12) # change to select by flight number...
    initial_ind = argmin(abs.((t.-initial_time))) # find the index of the initial time

    if ~isnothing(varg)
        # here we also are gonna need to check where things get inserted in case they are not in order...
        if !assume_monotonic # use data to figure out how and where to do insertions...
            # we need some way to get the local dimension from just a variable
        else
            vardata = combine_air_and_ground_data(vardata, vardatag, z_dim_num; insert_location=ground_indices) # append ground data with 0 as z bottom, loses labeling now though  (this puts a lot of faith im these 2 vars being the same size of having labels which we can't guarantee, no?)
        end
    end

    # select only the initial condition timestep, but keep that dimension around w/ []
    # note -- if not the initial condition, we should still only return initial condition to reference timestep no? (or i guess at least just from the initial condition to the end of the data we have...)
    if initial_condition # the timestep that is closest to the one we are supposed to force with (should be reference time - 12 hours)
        vardata  = selectdim(vardata, time_dim_num, [initial_ind])
        z_old    = selectdim(z_old, time_dim_num, [initial_ind])
    else # not init condition, so we'll truncate from init condition to end along time dimension
        vardata  = selectdim(vardata, time_dim_num, initial_ind:length(t_old))
        z_old    = selectdim(z_old, time_dim_num, initial_ind:length(t_old))
    end

    # reverse the z so it goes from ground to top) and matches the new grid we defined..
    z_old = reverse(z_old; dims=z_dim_num)
    vardata = reverse(vardata; dims=z_dim_num)

    # interpolate to new z
    vardata = var_to_new_coord(vardata, z_old,    z_dim_num; coord_new=z_new  , data=data)
    if initial_condition # no need to push further here since is init condition (maybe change later to return both?)
        return vardata
    end
    # create new time splines
    vardata = var_to_new_coord(vardata, t_old, time_dim_num; coord_new=nothing, data=data)

    return vardata
end



function drop_lat_lon(vardata; data=nothing, dims=nothing)
    # gotta do this at the end after we've made full use of our dimnames since our data is unlabeled
    # in general we could cheat in the future since the socrates order seems to always be lon lat lev time regardless of which dims exist...
    # i guess i also don't know what happens if we've turned the time dim into just a fcn -- in principle it's last so that shouldnt hurt

    if isnothing(dims)
        dims = tuple(collect(findfirst(x->x==dim,NC.dimnames(data["T"])) for dim in ["lat","lon"] )...) # base off temperature for now
    end

    vardata = dropdims(vardata, dims = dims)
    return vardata
end

"""
    insert_dims

add in new dimensions at position (never seemed to use , could get rid of this?)
"""
function insert_dims(data, ind; new_dim_sizes=[-1])
    sz_data = collect(size(data)) # array
    # insert!(sz_data,ind, new_dim_sizes[1]) # insert sz 1 at this location # can only do one dim
    splice!(sz_data,ind, [ new_dim_sizes..., sz_data[ind]]) # do this way cause splice itself deletes the current value so preserve it
    data =  reshape(data, sz_data...) # reshape
end

function calc_qg(Tg,pg; thermo_params)
    pvg           = TD.saturation_vapor_pressure.(thermo_params, Tg, TD.Liquid())
    molmass_ratio = TDP.molmass_ratio(thermo_params)
    qg            = (1 / molmass_ratio) .* pvg ./ (pg .- pvg) #Total water mixing ratio at surface , assuming saturation [ add source ]
    return qg
end


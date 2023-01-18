#= ... 
created jan 4 2023
@jbphyswx
... =#

##
# - Maybe swith to using DimensionalData or something so we can have derived quantities (i.e. ts <thermodynamic state>) as labeled arrays and not have to keep passing their parents around...
# -- seems to be a little cumbersome and no easy way to do netCDF reads and conversion to array/dataset like types
##

using .SOCRATES_Single_Column_Forcings: grid_heights # is there a way to automatically link this to the other file? can't use w/o importing the entire module locally to link things
 
# this feels like a lot to import and have in the project file, should i copy over and write my own parameters?
# can't even use the latest thermodynamics versino which i have local with that import...
# Pkg.develop(path="/home/jbenjami/Research_Schneider/CliMa/TurbulenceConvection.jl")
# @show(SOCRATES_Single_Column_Forcings)

# const TCP = SOCRATES_Single_Column_Forcings.TCP
param_set = SOCRATES_Single_Column_Forcings.Parameters.param_set
thermo_params = SOCRATES_Single_Column_Forcings.TCP.thermodynamics_params(param_set)

# thermo_params = param_set
# thermo_params = SOCRATES_Single_Column_Forcings.Parameters.thermodynamics_params(param_set)

# figure out these types... #

# default_new_z = vec(grid.zc.z) 
default_new_z = collect(0:100:0000) # for testing case 

#= Simple linear interpolation function, wrapping Dierckx (copied from TC.jl)  =#
function pyinterp(x, xp, fp)
    # @show(x,xp,fp)
    spl = Dierckx.Spline1D(xp, fp; k = 1)
    return spl(vec(x))
end


function swapaxes(a, dim1, dim2)
    """
    Swaps the two dimensions of an array.
    """
    perm = collect(1:ndims(a))
    perm[dim1] = dim2
    perm[dim2] = dim1
    return permutedims(a, perm)
end

function align_along_dimension(v,dim)
    """
    Assumes v has one non-singleton dimension which will be aligned along dim. if dim does not exist, v is expanded till dim exists
    """

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


function data_to_tsg(data; param_set=param_set)
    """
    Does this cause we don't have surface q values as of rn
    """
    thermo_params = SOCRATES_Single_Column_Forcings.TCP.thermodynamics_params(param_set)

    Tg = data["Tg"]
    pg = data["Ps"]

    # not sure what to do at surface, so assuming saturation at surface
    pvg           = TD.saturation_vapor_pressure.(thermo_params, Tg, TD.Liquid()) 
    molmass_ratio = TCP.molmass_ratio(param_set)

    qg = (1 / molmass_ratio) .* pvg ./ (pg .- pvg) #Total water mixing ratio at surface , assuming saturation

    tsg = TD.PhaseEquil_pTq.(thermo_params, pg, Tg, qg) 
    return tsg
end

function data_to_ts(data; do_combine_air_and_ground_data=false, param_set=param_set) # add type data is ncdataset
    """
    """
    thermo_params = SOCRATES_Single_Column_Forcings.TCP.thermodynamics_params(param_set)

    T = data["T"]
    q = data["q"]
    p = data["lev"]

    # put p/lev on correct dimension...
    p = align_along_dimension(p, get_dim_num("lev", data["T"]))

    ts = TD.PhaseEquil_pTq.(thermo_params, p, T, q)

    # @show(size(T),size(q),size(p),size(ts),"data_to_ts")

    if do_combine_air_and_ground_data
        concat_dim = get_dim_num("lev", T) # phaseequil reduces us down to lev... can't seem to just apply ufunc...
        tsg = data_to_tsg(data;param_set=param_set)
        ts = combine_air_and_ground_data(ts,tsg, concat_dim)
    end

    return ts
end

# convert pressure to altitude...
function lev_to_z( p::FT, T::FT, q::FT, pg::FT, Tg::FT, qg::FT; param_set=param_set, data=data) where {FT <: Real}
    """
    """
    thermo_params = TCP.thermodynamics_params(param_set) # can we replace this and only rely on Thermodynamics?
    ts            = TD.PhaseEquil_pTq.(thermo_params, p, T, q) 
    tsg           = TD.PhaseEquil_pTq.(thermo_params, pg, Tg, qtg) 
    return lev_to_z(ts, tsg; data=data)
end

function lev_to_z(ts, tsg; param_set=param_set, data=data )
    """
    ts is thermodynamic state
    tsg is thermodynamic state for ground
    """
    thermo_params = TCP.thermodynamics_params(param_set) # can we replace this and only rely on Thermodynamics?
    R_d       = TCP.R_d(param_set) # TD.Parameters.R_d(thermo_params)
    grav      = TCP.grav(param_set)  # TD.Parameters.grav(thermo_params)
    
    dimnames    = NC.dimnames(data["T"]) # use this as default cause calculating ts doesn't maintain dim labellings
    lev_dim_num = findfirst(x->x=="lev",dimnames)
    ldn         = lev_dim_num
    L           = size(ts,lev_dim_num)

    # @show(size(ts), size(tsg), ldn)
    # tsz       = cat(ts,tsg;dims=ldn)
    tsz       = combine_air_and_ground_data(ts,tsg,ldn;data=data,reshape_ground=true)

    Tvz       = TD.virtual_temperature.(thermo_params, tsz) # virtual temp, we havent returned these for now...
    pz        = TD.air_pressure.(thermo_params, tsz)
    Lz        = L+1 # cause we extended it using the ground...
    Tv_bar    = Statistics.mean((selectdim(Tvz, ldn, 1:Lz-1), selectdim(Tvz, lev_dim_num, 2:Lz)))
    p_frac    = selectdim( pz, lev_dim_num, 2:Lz) ./ selectdim( pz, ldn, 1:Lz-1)
    dz        = @. (R_d * Tv_bar / grav) * log(p_frac)
    #
    z         = reverse(cumsum(reverse(dz;dims=ldn);dims=ldn);dims=ldn) #cumsum(dz) # grid is already defined from  (from Grid.jl)
    s_sz      = collect(size(z)) # should now be same dims as T
    s_sz[ldn] = 1 # we just want to add a single slice of zeros for the ground -- these are all ocean cases so this should be fine for now...
    _z0       = zeros(Float64,s_sz...)
    z         = cat(z, _z0; dims=ldn)
    return z
end

function z_from_data(data; param_set=param_set)
    ts  = data_to_ts(data; param_set=param_set, do_combine_air_and_ground_data=false)
    tsg = data_to_tsg(data; param_set=param_set)
    return lev_to_z(ts, tsg; param_set=param_set, data=data)
end


function combine_air_and_ground_data(var,varg, concat_dim; data=nothing, reshape_ground=true)
    """
    var  : the data variable or string variable name for data in the air
    varg : the data variable or string variable name for data on the ground
    data : is a container from which the data can be accessed as a string in data[var(g)] form

    reshape_ground : expand the ground data to the same size as bottom slice of the air data (so you can pass in for example a single scalar or similar reduced dimension varg)
    """

    # if data is a string, read the data out from data (creates data `vardata` from `var` whether var is string or already is data)
    vardata    = isa(var ,String) ? data[var ] : var
    vardatag   = isa(varg,String) ? data[varg] : varg
    concat_dim = isa(concat_dim,String) ? get_dim_num(concat_dim, vardata) : concat_dim # string to numbers

    # @show(size(vardata),size(vardatag), concat_dim)

    # handle the ground data alignment
    if reshape_ground
        if isa(vardatag, Number) # allow for you to pass a single number and still concat
            sz_vardatag = collect(size(vardata)) # array
            sz_vardatag[concat_dim] = 1
            vardatag = fill(vardatag, sz_vardatag... ) # create array full with just this one value
        elseif isa(vardatag, NC.CFVariable) # vardatag should not have lev as a dimension so we need to add it in (and I guess check the others are in the same order)
            # @show(size(vardatag))
            dimnamesg = NC.dimnames(vardatag)
            dimnames  = NC.dimnames(vardata)
            # in same order as vardata just in case
            vardatag = reshape(vardatag, (size(vardatag)...,1)) # add trailing singleton for lev
            dimnamesg = [dimnamesg..., "lev"]
            vardatag = permutedims(vardatag, [findfirst(x->x==dim,dimnames) for dim in dimnamesg] ) # permute into the right order
            # @show(size(vardatag))
        elseif isa(vardatag, AbstractArray) # an unlabeled array, i think we can't guarantee then that the lev axis exists at all or what existing dimensions are so this just assumes everything is correct except the lev dimension being there
            # assume we need to add a new dimension at the same location as in the full array and order otherwise is preserved (quick check looks ok)
            sz_vardatag = collect(size(vardatag)) # array
            insert!(sz_vardatag,concat_dim,1) # insert sz 1 at this location 
            vardatag =  reshape(vardatag, sz_vardatag...) # reshape (just adds the singleton dimension in)
        end
    end

    # @show(size(vardata),size(vardatag), concat_dim)
    vardata = cat(vardata,vardatag;dims=concat_dim)
    return vardata
end


function get_dim_num(dim,nc_data=nothing)
    """
    get the dimension number of dim from nc_data
    """

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



function interp_along_dim(var, interp_dim, interp_dim_in; interp_dim_out=nothing,data=nothing, data_func = nothing, interp_dim_in_is_full_array=true, reshape_ground=true, verbose=false)
    " interpolation data 

    # var : the data to be interpolated
    # interp_dim: the name or dimension number along which we will do interpolation

    # interp_dim_in: the coordinate on which the data is currently aligned 
    # interp_dim_out: the coordinate on which we would like to interpolate the data to  

    the function will create splines along interp_dim_in...
    - if we set interp_dim_out, we actually evaluate the function along interp_dim at locations interp_dim_out,otherwise, we just return our unevaluated spline functions

    - data_func is a func to be applied to the raw data before it is processed, though perhaps it it most useful if interp_dim_out is unset and we are returning functions...
    - vectorize_in means your input is an array and you need to loop over it too (as opposed to just being a fixed template vector)

    To Do : decide types
    "

    # if data is a string, read the data out from data (creates data `vardata` from `var` whether var is string or already is data)
    vardata    = isa(var ,String) ? data[var ] : var


    # get the interpolation dimension and combine air and ground data
    interp_dim_num = get_dim_num(interp_dim, vardata)
    # vardata        = combine_air_and_ground_data(vardata ,vardatag, interp_dim_num; data=nothing, reshape_ground=true) # combine air and ground data... (also resolves strings to data)
 
    if !isnothing(data_func) # apply data_func if we need to
        vardata  = data_func(vardata)
    end

    # @show(interp_dim_in_is_full_array, (size(interp_dim_in),size(vardata)), (size(interp_dim_in)==size(vardata)))
 
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

function var_to_new_coord(var,coord_in, interp_dim; coord_new=nothing, data=nothing, data_func=nothing, kwargs...)
    """
    if coord_new is nothing, then will return functions...
    """
    vardata    = isa(var ,String) ? data[var ] : var
    if ~isnothing(data)
        interp_dim = isa(interp_dim, String) ? get_dim_num(interp_dim, data) : interp_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
    else
        interp_dim = isa(interp_dim, String) ? get_dim_num(interp_dim, vardata) : interp_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
    end

    # evaluate interp_dim_in_is_full_array based on the size of the input... interp_dim_in is full array false is much faster cause dont have to double loop in vectorization...
    return interp_along_dim( vardata, interp_dim, coord_in; interp_dim_out=coord_new, data=data, data_func=data_func, interp_dim_in_is_full_array=(size(coord_in)==size(vardata)), kwargs...)
end


function get_data_new_z_tfunc(var, varg, z_new, z_dim, time_dim; z_old=nothing, t_old=nothing, data=nothing, param_set=param_set, initial_condition = false)
    """
    Take data from our base setup, interpolate it to new z, then create time splines based on the t we have....
    to vectorize properly over z_new, it should be the same shape as vardata+vardata_g
    """

    # get the data and dimensions we're working on, 
    vardata    = isa(var ,String) ? data[var ] : var
    vardatag   = isa(varg,String) ? data[varg] : varg

    # combine air and ground data
    if ~isnothing(data)
        z_dim_num = isa(z_dim, String) ? get_dim_num(z_dim, vardata) : z_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
        time_dim_num = isa(time_dim, String) ? get_dim_num(time_dim, vardata) : time_dim
    else
        z_dim_num = isa(z_dim, String) ? get_dim_num(z_dim, vardata) : z_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
        time_dim_num = isa(time_dim, String) ? get_dim_num(time_dim, vardata) : time_dim
    end




    if isnothing(z_old)
        z_old = z_from_data(data; param_set=param_set); # uses ground value to create the old z, pads bottom w/ 0
    end
    if isnothing(t_old)
        t_old = data["tsec"][:] # check this unit was right in the files
    end

    vardata = combine_air_and_ground_data(vardata, vardatag, z_dim_num) # append ground data with 0 as z bottom, loses labeling now though  (this puts a lot of faith im these 2 vars being the same size of having labels which we can't guarantee, no?)

    # select only the first timestep, but keep that dimension around w/ []
    if initial_condition
        vardata  = selectdim(vardata, time_dim_num, [1])
        z_old    = selectdim(z_old, time_dim_num, [1])
    end

    # reverse the z so it goes from ground to top) and matches the new grid we defined..
    z_old = reverse(z_old; dims=z_dim_num)
    vardata = reverse(vardata; dims=z_dim_num)

    # interpolate to new z
    vardata = var_to_new_coord(vardata, z_old,    z_dim_num; coord_new=z_new  , data=data)
    if initial_condition # no need to push further here since is init condition
        return vardata
    end
    # create new time splines
    vardata_init = var_to_new_coord(vardata, t_old, time_dim_num; coord_new=nothing, data=data)

    return vardata
end



# =================== #


function main()
    data = load_data()
    nothing
 end
 

 if abspath(PROGRAM_FILE) == @__FILE__
     main()
 end


 function test()
    data = SOCRATES_Single_Column_Forcings.open_atlas_les_profile(9);
    new_z = data[:grid_data]
    data = data[:obs_data]


    T_new_z_one_shot = get_data_new_z_tfunc("T", "Tg", new_z,"lev","time", data=data, param_set=param_set)
    T_new_z_initial  = get_data_new_z_tfunc("T", "Tg", new_z,"lev","time", data=data, param_set=param_set,  initial_condition=true)

end

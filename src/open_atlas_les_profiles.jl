"""
    open_atlas_les_profile(flight_number::Int)

opens the files downloaded in download_atlas_les_profiles.jl
"""
function open_atlas_les_profile(flight_number::Int)
    FT = Float64 # idk how to pass on a type here without necessarily having to give a variable...
    atlas_dir = joinpath(dirname(@__DIR__), "Data", "Atlas_LES_Profiles") # doesn't seem to work to use @__DIR__ directly as a variable
    RF_num = "RF" * string(flight_number,pad=2)
    obs_filename  = joinpath(atlas_dir, RF_num * "_obs-based_SAM_input.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_obs-based_SAM_input.nc
    ERA5_filename = joinpath(atlas_dir, RF_num * "_ERA5-based_SAM_input_mar18_2022.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_ERA5-based_SAM_input_mar18_2022.nc
    grid_filename = joinpath(atlas_dir, RF_num * "_grd.txt")

    # @show(obs_filename, ERA5_filename, grid_filename)

    local obs_data, ERA5_data, grid_data # initialize cause try catch scope is closed
    # can't use do blocks here cause will close the files...
    try  # obs
        # obs_data = NC.Dataset(obs_filename,"r") do ds; ds; end
        obs_data = NC.Dataset(obs_filename,"r")
    catch e
        @warn e
        obs_data = nothing
    end

    try  # ERA5
        ERA5_data = NC.Dataset(ERA5_filename,"r")
    catch e
        @warn e
        ERA5_data = nothing
    end

    try  # grid
        grid_data = vec(readdlm(grid_filename, FT)) # is a txt file...
    catch e
        @warn e
        grid_data = nothing
    end

    return (; obs_data = obs_data, ERA5_data = ERA5_data, grid_data = grid_data)

end
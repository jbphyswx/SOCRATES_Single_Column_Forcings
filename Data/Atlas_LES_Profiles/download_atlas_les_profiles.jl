#=
Download LES forcing data from https://atmos.uw.edu/~ratlas/SOCRATES-LES-cases.html
=#

import NCDatasets as NC
import SOCRATESSingleColumnForcings as SSCF

cases = collect(keys(SSCF.grid_heights))

thisdir = @__DIR__ # doesn't seem to work to use @__DIR__ directly as a variable

# download forcings for each flight

"""
    download_atlas_les_inputs(;cases=cases)

"""
function download_atlas_les_inputs(;cases=cases)
    for flight in cases # the socrates flight numbers
        RF_num = "RF" * string(flight,pad=2)
        obs_filename  = RF_num * "_obs-based_SAM_input.nc" # e.g. https://atmos.uw.edu/~ratlas/RF12_obs-based_SAM_input.nc
        ERA5_filename = RF_num * "_ERA5-based_SAM_input_mar18_2022.nc" # e.g. https://atmos.uw.edu/~ratlas/RF12_ERA5-based_SAM_input_mar18_2022.nc

        # download obs forcing
        try
            obs_savepath = joinpath(thisdir, "Input_Data", obs_filename)
            download("https://atmos.uw.edu/~ratlas/"*obs_filename, obs_savepath)
            @info "Found $("https://atmos.uw.edu/~ratlas/"*obs_filename)"
        catch e
            @warn "Did not find $("https://atmos.uw.edu/~ratlas/"*obs_filename)"
            # @warn e
        end

        # download ERA5 forcing
        try
            ERA5_savepath = joinpath(thisdir, "Input_Data", ERA5_filename)
            download("https://atmos.uw.edu/~ratlas/"*ERA5_filename, ERA5_savepath)
            @info "Found $("https://atmos.uw.edu/~ratlas/"*ERA5_filename)"
        catch e
            @warn "Did not find $("https://atmos.uw.edu/~ratlas/"*ERA5_filename)"
            # @warn e
        end

        #download grid file (is same for both era and obs forcings)
        grid_height = SSCF.grid_heights[flight]
        try
            download("https://atmos.uw.edu/~ratlas/"*string(grid_height)*"level-grd.txt",  thisdir*"/"*RF_num*"_grd.txt")
            @info "Found https://atmos.uw.edu/~ratlas/"*string(grid_height)*"level-grd.txt"
        catch e
            @warn "Did not find https://atmos.uw.edu/~ratlas/"*string(grid_height)*"level-grd.txt"
        end
    end
end





"""
    download_atlas_les_profiles(;cases=cases) 

    These all have the same fiilename so we have to manually append the flight number to the filename

"""
function download_atlas_les_outputs(;cases=cases)
    for flight in cases # the socrates flight numbers
        RF_num = "RF" * string(flight,pad=2)
        
        obs_filename  = RF_num * "_output/obs/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc"  # e.g. https://atmos.uw.edu/~ratlas/RF13_output/obs/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc
        ERA5_filename = RF_num * "_output/era5/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" # e.g. https://atmos.uw.edu/~ratlas/RF12_output/era5/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc

        # download obs model output
        try
            obs_savepath = joinpath(thisdir,"Output_Data", RF_num*"_Obs_" * split(obs_filename,"/")[end]) # just the filename, not any paths
            download("https://atmos.uw.edu/~ratlas/"*obs_filename, obs_savepath)
            @info "Found $("https://atmos.uw.edu/~ratlas/"*obs_filename)"
        catch e
            @warn "Did not find $("https://atmos.uw.edu/~ratlas/"*obs_filename)"
            # @warn e
        end

        # download ERA5 forcing
        try
            ERA5_savepath = joinpath(thisdir, "Output_Data", RF_num*"_ERA5_" * split(ERA5_filename,"/")[end]) # just the filename, not any paths
            download("https://atmos.uw.edu/~ratlas/"*ERA5_filename, ERA5_savepath)
            @info "Found $("https://atmos.uw.edu/~ratlas/"*ERA5_filename)"
        catch e
            @warn "Did not find $("https://atmos.uw.edu/~ratlas/"*ERA5_filename)"
            # @warn e
        end
    end
end

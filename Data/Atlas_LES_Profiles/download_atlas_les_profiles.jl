#=
Download LES forcing data from https://atmos.uw.edu/~ratlas/SOCRATES-LES-cases.html
=#

import NCDatasets as NC
import SOCRATESSingleColumnForcings as SSCF

cases = collect(keys(SSCF.grid_heights))

thisdir = @__DIR__ # doesn't seem to work to use @__DIR__ directly as a variable

# download forcings for each flight

"""
    download_atlas_les_profiles(;cases=cases)

"""
function download_atlas_les_profiles(;cases=cases)
    for flight in cases # the socrates flight numbers
        RF_num = "RF" * string(flight,pad=2)
        obs_filename  = RF_num * "_obs-based_SAM_input.nc" # e.g. https://atmos.uw.edu/~ratlas/RF12_obs-based_SAM_input.nc
        ERA5_filename = RF_num * "_ERA5-based_SAM_input_mar18_2022.nc" # e.g. https://atmos.uw.edu/~ratlas/RF12_ERA5-based_SAM_input_mar18_2022.nc

        # download obs forcing
        try
            obs_savepath = thisdir*"/"*obs_filename
            download("https://atmos.uw.edu/~ratlas/"*obs_filename, obs_savepath)
            @warn "Found $("https://atmos.uw.edu/~ratlas/"*obs_filename)"

        catch e
            @warn "Did not find $("https://atmos.uw.edu/~ratlas/"*obs_filename)"
            # @warn e
        end

        # download ERA5 forcing
        try
            ERA5_savepath = thisdir*"/"*ERA5_filename
            download("https://atmos.uw.edu/~ratlas/"*ERA5_filename, ERA5_savepath)
            @warn "Found $("https://atmos.uw.edu/~ratlas/"*ERA5_filename)"
        catch e
            @warn "Did not find $("https://atmos.uw.edu/~ratlas/"*ERA5_filename)"
            # @warn e
        end

        #download grid file (is same for both era and obs forcings)
        grid_height = SSCF.grid_heights[flight]
        try
            download("https://atmos.uw.edu/~ratlas/"*string(grid_height)*"level-grd.txt",  thisdir*"/"*RF_num*"_grd.txt")
            @warn "Found https://atmos.uw.edu/~ratlas/"*string(grid_height)*"level-grd.txt"
        catch e
            @warn "Did not find https://atmos.uw.edu/~ratlas/"*string(grid_height)*"level-grd.txt"
        end
    end
end

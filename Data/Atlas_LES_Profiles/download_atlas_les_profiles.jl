"""
Download LES forcing data from https://atmos.uw.edu/~ratlas/SOCRATES-LES-cases.html
"""

export grid_heights
export cases

PKGDIR="/home/jbenjami/Research_Schneider/CliMa/SOCRATES_Single_Column_Forcings.jl";
empty!(DEPOT_PATH); push!(DEPOT_PATH,PKGDIR*"/.julia_depot"); 
import NCDatasets as NC

# Two cases with shallow cloud-topped boundary layers, RF12 and RF13, are run on a 192-level vertical grid. 
# The other four cases have clouds extending through deeper boundary layers; they are run on a 320-level vertical grid.
grid_heights = Dict(
     1 => 320,
     9 => 320,
    10 => 320,
    11 => 320,
    12 => 192, 
    13 => 192,
)

cases = collect(keys(grid_heights))

thisdir = @__DIR__ # doesn't seem to work to use @__DIR__ directly as a variable

# download forcings for each flight

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
        grid_height = grid_heights[flight]
        try
            download("https://atmos.uw.edu/~ratlas/"*string(grid_height)*"level-grd.txt",  thisdir*"/"*RF_num*"_grd.txt")
            @warn "Found https://atmos.uw.edu/~ratlas/"*string(grid_height)*"level-grd.txt"
        catch e
            @warn "Did not find https://atmos.uw.edu/~ratlas/"*string(grid_height)*"level-grd.txt"
        end
    end
end

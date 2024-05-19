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
function download_atlas_les_inputs(; cases = cases)
    for flight in cases # the socrates flight numbers
        RF_num = "RF" * string(flight, pad = 2)
        obs_filename = RF_num * "_obs-based_SAM_input.nc" # e.g. https://atmos.uw.edu/~ratlas/RF12_obs-based_SAM_input.nc
        ERA5_filename = RF_num * "_ERA5-based_SAM_input_mar18_2022.nc" # e.g. https://atmos.uw.edu/~ratlas/RF12_ERA5-based_SAM_input_mar18_2022.nc

        # download obs forcing
        try
            obs_savepath = joinpath(thisdir, "Input_Data", obs_filename)
            download("https://atmos.uw.edu/~ratlas/" * obs_filename, obs_savepath)
            @info "Found $("https://atmos.uw.edu/~ratlas/"*obs_filename)"
        catch e
            @warn "Did not find $("https://atmos.uw.edu/~ratlas/"*obs_filename)"
            # @warn e
        end

        # download ERA5 forcing
        try
            ERA5_savepath = joinpath(thisdir, "Input_Data", ERA5_filename)
            download("https://atmos.uw.edu/~ratlas/" * ERA5_filename, ERA5_savepath)
            @info "Found $("https://atmos.uw.edu/~ratlas/"*ERA5_filename)"
        catch e
            @warn "Did not find $("https://atmos.uw.edu/~ratlas/"*ERA5_filename)"
            # @warn e
        end

        #download grid file (is same for both era and obs forcings)
        grid_height = SSCF.grid_heights[flight]
        try
            download(
                "https://atmos.uw.edu/~ratlas/" * string(grid_height) * "level-grd.txt",
                thisdir * "/" * RF_num * "_grd.txt",
            )
            @info "Found https://atmos.uw.edu/~ratlas/" * string(grid_height) * "level-grd.txt"
        catch e
            @warn "Did not find https://atmos.uw.edu/~ratlas/" * string(grid_height) * "level-grd.txt"
        end
    end
end





SOCRATES_flight_observations_Box_links = Dict( # The raw observational data from the SOCRATES flights, for ready access on any machine...
    # They're big files so not storing them in the git repo...
    1 => "https://caltech.box.com/shared/static/4sh03gadjq6acary69qamgjbmynfio2f.nc",
    2 => "https://caltech.box.com/shared/static/urhtegy7dccnp8hav4my90kzrr3e4wmt.nc",
    3 => "https://caltech.box.com/shared/static/xp4b4p2ef523bqmzuejbecmruu8daj6s.nc",
    4 => "https://caltech.box.com/shared/static/f41qh4jjlnvs08y676hhq60epcta2oip.nc",
    5 => "https://caltech.box.com/shared/static/ty0h3mxrj6myyun3qstfoj1ef0o8009m.nc",
    6 => "https://caltech.box.com/shared/static/e3yk43xfwyxyix0yqsh53utpmjynzkzr.nc",
    7 => "https://caltech.box.com/shared/static/34bz9z9hlolqwmer02g2qz4uu13qn13s.nc",
    8 => "https://caltech.box.com/shared/static/pua0gf97772xn5cdq256bztriu4q5kcz.nc",
    9 => "https://caltech.box.com/shared/static/1vdwhg0oiwtzd21vqwzt5sm6zsh8yng9.nc",
    10 => "https://caltech.box.com/shared/static/cw40hwxoegcmjxksja4gxqfpa22mwj9p.nc",
    11 => "https://caltech.box.com/shared/static/cz7qhv4ub5kbsxloxydf2kyykd9e5xdn.nc",
    12 => "https://caltech.box.com/shared/static/mugnp8iqfcccksxd8uojqbufcn8waexh.nc",
    13 => "https://caltech.box.com/shared/static/qzj14slfrknboi4iflen8t5gqfgmi595.nc",
    14 => "https://caltech.box.com/shared/static/kxcaogfc1m5eoopb5z3gedffv8gjizax.nc",
    15 => "https://caltech.box.com/shared/static/pfp8egkqcxhbyi1bxqfgk28g7neih7d8.nc",
)



"""
    download_atlas_les_profiles(;cases=cases) 

    These all have the same fiilename so we have to manually append the flight number to the filename

"""
function download_atlas_les_outputs(; cases = cases)
    for flight in cases # the socrates flight numbers
        RF_num = "RF" * string(flight, pad = 2)

        obs_filename = RF_num * "_output/obs/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc"  # e.g. https://atmos.uw.edu/~ratlas/RF13_output/obs/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc
        ERA5_filename = RF_num * "_output/era5/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" # e.g. https://atmos.uw.edu/~ratlas/RF12_output/era5/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc

        # download obs model output
        try
            obs_savepath = joinpath(thisdir, "Output_Data", RF_num * "_Obs_" * split(obs_filename, "/")[end]) # just the filename, not any paths
            download("https://atmos.uw.edu/~ratlas/" * obs_filename, obs_savepath)
            @info "Found $("https://atmos.uw.edu/~ratlas/"*obs_filename)"
        catch e
            @warn "Did not find $("https://atmos.uw.edu/~ratlas/"*obs_filename)"
            # @warn e
        end

        # download ERA5 forcing
        try
            ERA5_savepath = joinpath(thisdir, "Output_Data", RF_num * "_ERA5_" * split(ERA5_filename, "/")[end]) # just the filename, not any paths
            download("https://atmos.uw.edu/~ratlas/" * ERA5_filename, ERA5_savepath)
            @info "Found $("https://atmos.uw.edu/~ratlas/"*ERA5_filename)"
        catch e
            @warn "Did not find $("https://atmos.uw.edu/~ratlas/"*ERA5_filename)"
            # @warn e
        end
    end
end

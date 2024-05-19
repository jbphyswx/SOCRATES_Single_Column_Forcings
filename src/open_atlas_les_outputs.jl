
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
open_atlas_les_output(flight_number::Int)

opens the files downloaded in download_atlas_les_profiles.jl
"""
function open_atlas_les_output(flight_number::Int)
    FT = Float64 # idk how to pass on a type here without necessarily having to give a variable...
    atlas_dir = joinpath(dirname(@__DIR__), "Data", "Atlas_LES_Profiles", "Output_Data") # doesn't seem to work to use @__DIR__ directly as a variable
    RF_num = "RF" * string(flight_number,pad=2)
    obs_filename  = joinpath(atlas_dir, RF_num * "_Obs_" *  "SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_obs-based_SAM_input.nc
    ERA5_filename = joinpath(atlas_dir, RF_num * "_ERA5_" *  "SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_ERA5-based_SAM_input_mar18_2022.nc
    grid_filename = joinpath(atlas_dir, RF_num * "_grd.txt")

    # @show(obs_filename, ERA5_filename, grid_filename)

    # check if file exists and if not, download it
    if !isfile(obs_filename)
        # download(SOCRATES_flight_observations_Box_links[flight_number], obs_filename) # box link
        download_atlas_les_outputs(;cases=[flight_number]) # Atlas website
    end
    if !isfile(ERA5_filename)
        # download(SOCRATES_flight_observations_Box_links[flight_number], ERA5_filename) # box link
        download_atlas_les_outputs(;cases=[flight_number]) # Atlas website
    end


   # The data is quite large so we'll try to load it from Box

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


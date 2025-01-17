
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
open_atlas_les_output(flight_number::Int, forcing_type::Symbol; open_files::Bool = true, include_grid::Bool = true)

opens the files downloaded in download_atlas_les_profiles.jl
"""
function open_atlas_les_output(
    flight_number::Int,
    forcing_type::Symbol;
    open_files::Bool = true,
    include_grid::Bool = true,
)
    FT = Float64 # idk how to pass on a type here without necessarily having to give a variable...
    atlas_dir = joinpath(dirname(@__DIR__), "Data", "Atlas_LES_Profiles", "Output_Data") # doesn't seem to work to use @__DIR__ directly as a variable
    RF_num = "RF" * string(flight_number, pad = 2)

    if forcing_type == :obs_data
        data_filename = joinpath(atlas_dir, RF_num * "_Obs_" * "SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_obs-based_SAM_input.nc
    elseif forcing_type == :ERA5_data
        data_filename = joinpath(atlas_dir, RF_num * "_ERA5_" * "SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_ERA5-based_SAM_input_mar18_2022.nc
    else
        error("forcing_type must be :obs_data or :ERA5_data")
    end


    # check if file exists and if not, download it
    # The data is quite large so we'll try to load it from Box
    if !isfile(data_filename)
        download_atlas_les_outputs(; cases = [flight_number], forcing_type = forcing_type) # Atlas website
    end


    data = isfile(data_filename) ? (open_files ? NC.Dataset(data_filename, "r") : data_filename) : nothing

    if include_grid
        grid_filename = joinpath(atlas_dir, RF_num * "_grd.txt")
        grid_data = isfile(grid_filename) ? (open_files ? vec(readdlm(grid_filename, FT)) : grid_filename) : nothing
        return NamedTuple{(forcing_type, :grid_data)}((data, grid_data))
    else
        return NamedTuple{(forcing_type,)}((data,))
    end
end


function open_atlas_les_output(flight_number::Int; open_files::Bool = true, include_grid::Bool = true)
    FT = Float64 # idk how to pass on a type here without necessarily having to give a variable...
    atlas_dir = joinpath(dirname(@__DIR__), "Data", "Atlas_LES_Profiles", "Output_Data") # doesn't seem to work to use @__DIR__ directly as a variable
    RF_num = "RF" * string(flight_number, pad = 2)
    obs_filename = joinpath(atlas_dir, RF_num * "_Obs_" * "SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_obs-based_SAM_input.nc
    ERA5_filename = joinpath(atlas_dir, RF_num * "_ERA5_" * "SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_ERA5-based_SAM_input_mar18_2022.nc

    # check if file exists and if not, download it
    # The data is quite large so we'll try to load it from Box
    if !isfile(obs_filename)
        download_atlas_les_outputs(; cases = [flight_number], forcing_type = :obs_data) # Atlas website
    end
    if !isfile(ERA5_filename)
        download_atlas_les_outputs(; cases = [flight_number], forcing_type = :ERA5_data) # Atlas website
    end

    obs_data = isfile(obs_filename) ? (open_files ? NC.Dataset(obs_filename, "r") : obs_filename) : nothing
    ERA5_data = isfile(ERA5_filename) ? (open_files ? NC.Dataset(ERA5_filename, "r") : ERA5_filename) : nothing

    if include_grid
        grid_filename = joinpath(atlas_dir, RF_num * "_grd.txt")
        grid_data = isfile(grid_filename) ? (open_files ? vec(readdlm(grid_filename, FT)) : grid_filename) : nothing
        return (; obs_data, ERA5_data, grid_data)
    else
        return (; obs_data, ERA5_data)
    end
end

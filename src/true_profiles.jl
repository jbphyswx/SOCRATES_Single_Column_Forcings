"""
This is a program to generate true profiles of variables from the data...

NSF/NCAR GV HIAPER 2D-S Particle Size Distribution (PSD) Product Data
https://data.eol.ucar.edu/dataset/552.047

# storing atlas processed files at https://caltech.app.box.com/folder/188886449133

"""

using Downloads # move to main

# store here cause they're too big... files are in files at https://caltech.app.box.com/folder/188886449133, could give em nicer names and stuff but would need to figure out an API for downloading...

# These don't exactly work and I don't wanna fight with the Box or Drive API -- could use Git LFS but does that mean every client needs it installed? I guess ask what they do for other stuff... on sampo whre my git is local could at least use https://anaconda.org/conda-forge/git-lfss
# (use box shared static links! -- see link settings on box link)
processed_Atlas_flight_data = Dict(
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

# These don't 

function open_processed_Atlas_flight_data(flight_number::Int; file_save_location=nothing, cleanup_file::Bool=true)
"""
"""

# file_save_location = isnothing(file_save_location) ? "/tmp/"*string(uuid1())*".nc"  : file_save_location
file_save_location = isnothing(file_save_location) ? tempname()*".nc"  : file_save_location

Downloads.download(processed_Atlas_flight_data[flight_number], file_save_location)

data = NCDatasets.Dataset(file_save_location, "r")

if cleanup_file # note /tmp would get cleaned anywaybut say you saved somewhere else, would need to load into memory though
end

return data

end


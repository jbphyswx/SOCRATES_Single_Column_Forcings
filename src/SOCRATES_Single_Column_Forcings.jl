module SOCRATES_Single_Column_Forcings

# put these first so are avail everywhere (removed cause doesn't work if you're say loading the package from github)
# THIS_DIR="/home/jbenjami/Research_Schneider/CliMa/SOCRATES_Single_Column_Forcings.jl";
# empty!(DEPOT_PATH); push!(DEPOT_PATH,THIS_DIR*"/.julia_depot"); 
using Pkg
Pkg.Registry.add(Pkg.RegistrySpec(url="https://github.com/jbphyswx/MyRegistry")) # see https://discourse.julialang.org/t/more-informative-project-toml-or-partial-manifest-toml/92804/5 

# Pkg.develop(path="/home/jbenjami/Research_Schneider/CliMa/Thermodynamics.jl") # do i still need this?
# Pkg.add(url="https://github.com/CliMA/Thermodynamics.jl#jb/non_eq_moisture") # this still worked to compile this packge, but somethign goes wrong when using it in TC.jl
# Pkg.add(name="Thermodynamics",rev="jb/non_eq_moisture") # doesn't seem your localregistry specifies branch so i guess you'd have to do this... (i think having this uncommented makes an infinite recursion when adding from the environment  )
# Pkg.add(url="https://github.com/CliMA/Thermodynamics.jl.git",rev="jb/non_eq_moisture") # force track git branch , no registry name clash i hope, see https://github.com/JuliaLang/Pkg.jl/issues/681#issuecomment-415190878 (seems to make precompilation hang and stall..... idk why -- maybe cause i removed from localregistry?)
# Pkg.develop(url="https://github.com/CliMA/Thermodynamics.jl.git#jb/non_eq_moisture")  # this does work but adds another folder...
# i think you just need to make sure the manifest has the tracked branch and just don't do all that here...

import Thermodynamics as TD
import NCDatasets as NC
using DelimitedFiles
using Statistics
using Dierckx
# import .Parameters as TCP


# include our files
include("Parameters.jl")
using .Parameters
const TCP = SOCRATES_Single_Column_Forcings.Parameters # import doesnt seem to work (has to go first to expose this to the other files)
FT = Float64
flight_numbers = [1, 9, 10, 11, 12 ,13]
forcing_types  = [:obs_data, :ERA5_data] # maybe change these to [:obs,:ERA5] later? would need to mirror in Cases.jl in TC.jl

include("../Data/Atlas_LES_Profiles/download_atlas_les_profiles.jl")
include("../Data/Atlas_LES_Profiles/open_atlas_les_profiles.jl")
include("helper_functions.jl")
include("process_case.jl")

end
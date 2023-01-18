module SOCRATES_Single_Column_Forcings

# put these first so are avail everywhere
THIS_DIR="/home/jbenjami/Research_Schneider/CliMa/SOCRATES_Single_Column_Forcings.jl";
empty!(DEPOT_PATH); push!(DEPOT_PATH,THIS_DIR*"/.julia_depot"); 
using Pkg
Pkg.develop(path="/home/jbenjami/Research_Schneider/CliMa/Thermodynamics.jl")
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

include("../Data/Atlas_LES_Profiles/download_atlas_les_profiles.jl")
include("../Data/Atlas_LES_Profiles/open_atlas_les_profiles.jl")
include("helper_functions.jl")



end
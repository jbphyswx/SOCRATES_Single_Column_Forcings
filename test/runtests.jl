using Test
import SOCRATESSingleColumnForcings as SSCF
import CLIMAParameters as CP
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

@testset "SOCRATESSingleColumnForcings" begin
    FT = Float64
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias");
    aliases = string.(fieldnames(TDP.ThermodynamicsParameters));
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics");
    thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...);

    data = SSCF.open_atlas_les_profile(9);
    new_z = data[:grid_data]
    data = data[:obs_data]
    # dTdt_hadv  = get_data_new_z_t(dTdt_hadv , new_z, z_dim_num,time_dim_num; z_old = z_old[:ERA5_data], data=data[:ERA5_data], thermo_params,  initial_condition)
    # var, z_new, z_dim, time_dim
    # T_new_z_t     = SSCF.get_data_new_z_t("T", new_z,"lev","time"; data, thermo_params)
    # T_new_z_init  = SSCF.get_data_new_z_t("T", new_z,"lev","time"; data, thermo_params,  initial_condition=true)
end

using Test
import SOCRATESSingleColumnForcings as SSCF
import CLIMAParameters as CP # use CLIMAParameters = "0.7, 0.8, 0.9, 0.10"
# import ClimaParams as CPP # would using this trouble w/ TC.jl? it's a different uuid technically..., use ClimaParams = "0.10"
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

@testset "SOCRATESSingleColumnForcings" begin
    FT = Float64

    toml_dict = CP.create_toml_dict(FT; dict_type = "alias") # CP 0.7 and below, Thermodynamics 0.11 and above
    aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)

    # toml_dict = CP.create_toml_dict(FT;) # CP 0.8 and up need to figure out
    # thermo_params = TDP.ThermodynamicsParameters(toml_dict)

    # toml_dict = CPP.create_toml_dict(FT;) # ClimaParams 0.10+ (for use w/ cloudmicrophysics 0.18+)
    # thermo_params = TDP.ThermodynamicsParameters(toml_dict)

    data = SSCF.open_atlas_les_input(9)
    new_z = data[:grid_data]
    data = data[:obs_data]
    # dTdt_hadv  = get_data_new_z_t(dTdt_hadv , new_z, z_dim_num,time_dim_num; z_old = z_old[:ERA5_data], data=data[:ERA5_data], thermo_params,  initial_condition)
    # var, z_new, z_dim, time_dim
    # T_new_z_t     = SSCF.get_data_new_z_t("T", new_z,"lev","time"; data, thermo_params)
    # T_new_z_init  = SSCF.get_data_new_z_t("T", new_z,"lev","time"; data, thermo_params,  initial_condition=true)

    # these no longer work w/o output_data downloaded
    data = SSCF.process_case(9; thermo_params = thermo_params, use_LES_output_for_z = false)
    @show(data)

    data = SSCF.process_case(9; thermo_params = thermo_params, initial_condition = true, use_LES_output_for_z = false)
    @show(data)

    data =
        SSCF.process_case(9; thermo_params = thermo_params, surface = "reference_state", use_LES_output_for_z = false)
    @show(data)

    data = SSCF.process_case(
        9;
        thermo_params = thermo_params,
        surface = "surface_conditions",
        use_LES_output_for_z = false,
    )
    @show(data)

end

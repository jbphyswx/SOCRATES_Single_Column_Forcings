# SOCRATES_Single_Column_Forcings

This repository includes scripts to take in SOCRATES flight data and create both 2D [lev,time] datasets and 1D [lev] columns of time-interpoloation functions to force the CliMa EDMF Turbulence-Convection single-column model.

The profiles are created off the data used to force LES simulations by Atlas, 2020 (https://doi.org/10.1029/2020MS002205), some of her scripts are included in /Rachel_Atlas_Scripts.

The output data are too large to store here and so will be stored in box. They 1D functions are interpolated to the grid that Atlas used, but functions are provided to interpolate to your own grid/time and return either raw data in time or a column of spline functions...

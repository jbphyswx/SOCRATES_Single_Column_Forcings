import numpy as np
from netCDF4 import Dataset, num2date
import glob
import subset_flight_table as sft
from scipy import interpolate, stats
import pandas
from datetime import date, datetime, timedelta
from subprocess import call

Mair=28.97
R=8.314

dp=2 #Choose pressure bin width in millibars
p_bins=np.arange(400,1010,dp)
p_bins_mid=p_bins[0:-1]+np.diff(p_bins)/2

#Replace with the flights that you are interested in
for flight_no in [1,9,10,12,13]:

    if flight_no<10:
      flight='rf0'+str(flight_no)
    else:
      flight='rf'+str(flight_no)
  
    #Create a new directory to store binned data- replace the path!
    call("rm -r ~/for_Jordan/"+\
         flight+'_'+str(dp)+'mb',\
         shell=True)
  
    call("mkdir ~/for_Jordan/"+\
         flight+'_'+str(dp)+'mb',\
         shell=True)
  
    #Replace with path to your SOCRATES files
    path = '~/for_Jordan/'
    if flight_no<10:
      file_search=glob.glob(path+'RF0'+str(flight_no)+'*may23.nc')
    else:
      file_search=glob.glob(path+'RF'+str(flight_no)+'*may23.nc')
 
    #Load dataset- can replace with xarray
    low_rate_filename=file_search[0]
    low_rate=Dataset(low_rate_filename)
  
    alt=low_rate.variables['GGALT'][:]
    pres=low_rate.variables['PSXC'][:]
    T=low_rate.variables['ATX'][:]
    theta=low_rate.variables['THETA'][:]
    #compute density to convert LWC
    rho=pres*1e-1*Mair/((T+273.15)*R)
  
    #Load variables that will be binned
    T=low_rate.variables['ATX'][:]
    w_obs=low_rate.variables['WIC'][:]
    qv=low_rate.variables['MR'][:]

    #The cloud flag is something that I added to the files
    #There is a description of this flag in the netCDF files
    cloud_flag=low_rate.variables['CLOUD_FLAG'][:]
    cloud_flag[cloud_flag==2]=0
  
    #Load liquid/condensed water content from King probe/CDP 
    lw_king=low_rate.variables['PLWCC'][:]
    cw_cdp=low_rate.variables['PLWCD_RWIO'][:]
  
    #Find start and end times for all slantwise profiles in the flight
    #Consult subset_flight_table.py for more info
    subset=sft.subset_flight_table(['U','D'],[flight],120,\
           MakeArray=True,\
           NoSawtooth=True,modules=True)
    starts=subset[:,4].astype(int)
    ends=subset[:,5].astype(int)
    
    #Mask bad values of qv
    qv=np.ma.masked_where(qv<-1000.,qv)
  
    profile_number=0
  
    #Loop through flight legs and subset and bin variables
    for w in range(len(subset)):
 
        #Make sure that the end of the flight leg is not later
        # than the end of the file
        if ends[w]>len(T):
          ends[w]=len(T)

        #subset pressure  
        pres_sub=pres[starts[w]:ends[w]]

        #subset and bin altitude
        alt_sub=alt[starts[w]:ends[w]]
        [z_binned,_,_]=stats.binned_statistic(pres_sub,\
          alt_sub,statistic='median',bins=p_bins)
  
        #subset and bin temperature
        T_sub=np.squeeze(T[starts[w]:ends[w]])
        [T_binned,_,_]=stats.binned_statistic(pres_sub,\
          T_sub,statistic='median',bins=p_bins)
  
        #subset and bin temperature
        theta_sub=np.squeeze(theta[starts[w]:ends[w]])
        [theta_binned,_,_]=stats.binned_statistic(pres_sub,\
          theta_sub,statistic='median',bins=p_bins)
  
        #subset and bin humidity
        qv_sub=qv[starts[w]:ends[w]]
  
        #Make sure not to bin NaNs
        try:
          qv_index=np.argwhere(np.isnan(qv_sub)==False)
          [qv_binned,_,_]=stats.binned_statistic(np.squeeze\
                          (pres_sub[qv_index]),\
                          np.squeeze(qv_sub[qv_index]),\
                          statistic='median',bins=p_bins)
        except:
          qv_binned=0.0*p_bins_mid
  
        #Use cloud flag to only use in-cloud sampled for binning LWC
        cloud_flag_sub=cloud_flag[starts[w]:ends[w]]
        in_cloud=np.squeeze(np.where(cloud_flag_sub==1))+starts[w]

        pres_in_cloud=np.squeeze(pres[in_cloud])  
        king_sub=np.squeeze(lw_king[in_cloud])
        cdp_sub=np.squeeze(cw_cdp[in_cloud])
        rho_sub=np.squeeze(rho[in_cloud])

        #Make a flag for where LWC data is neither zero nor missing
        good_cdp=np.nonzero(cdp_sub)
        good_king=np.nonzero(king_sub)
 
        #Use CDP unless it wasn't working, then use King Probe
        #If there are no nonzero values for either instrument, 
        #all data is clear sky and lwc is set to NaN

        if np.shape(good_cdp)[1]>1: #CDP is working and there is cloud
          [water_binned,_,_]=stats.binned_statistic(pres_in_cloud\
            [good_cdp],cdp_sub[good_cdp]/rho[good_cdp],\
            statistic='median',bins=p_bins)
        elif np.shape(good_king)[1]>1: #CDP was not working and there is cloud
          [water_binned,_,_]=stats.binned_statistic(pres_in_cloud\
            [good_king],king_sub[good_king]/rho[good_king],\
            statistic='median',bins=p_bins)
        else: #There is no cloud
          water_binned=np.zeros(len(p_bins)-1)+np.nan 
 

	#Save binned profiles to text files
        np.savetxt('~/for_Jordan/'+flight+'_'+str(dp)+'mb/'+flight+'_'+\
                   str(profile_number)+\
                   '.txt',np.column_stack([T_binned,qv_binned,\
                   water_binned,z_binned,theta_binned]))
        profile_number=profile_number+1
    
    #Close file    
    low_rate.close() 

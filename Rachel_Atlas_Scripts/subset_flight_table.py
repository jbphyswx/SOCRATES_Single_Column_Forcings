"""
Created Apr 2018

@author: cloud atlas
"""

#This function searches for flight legs according to the leg type, 
#research flight and a minimum leg duration.  If MakeFile=True, the
#code will make a text file with the flight leg data.  
#If MakeArray=True, the function will return the leg data as an array.
#The function will always print the flight leg data to the terminal.

#Example function calls are:
#subset_flight_table(['A'],['rf01'],60)
#subset_flight_table(['A','B'],['rf01','rf02'],60) 
#subset_flight_table(['A','B'],['all'],60)
#subset_flight_table(['A'],['all',],60,MakeFile=True) 
#subset_flight_table(['A'],['all'],60,MakeFile=True,MakeArray=True)

def subset_flight_table(leg_types,flights,min_time,MakeFile=False,MakeArray=False):

  import numpy as np
  import pandas as pd
  import os

  #DO NOT CHANGE THIS PATH- the most updated leg table will be here 
  path = '~/for_jordan/'

  #T=Tarmac, F=Ferry, U=Up sounding, D=Down sounding, C=In-cloud level   #leg, B=Below cloud level leg, A=Above cloud level leg, S=Sawtooth,
  #CL=Clear sky
  leg_type_options=['T','F','U','D','C','B','A','S','CL']
 
  #Specify an array of any number of flight legs or 'all'
  flight_options=['rf01','rf02','rf03','rf04','rf05','rf06','rf07','rf08','rf09','rf10','rf11','rf12','rf13','rf14','rf15','all']

  #Check that all requested leg types and flights exist
  n=0
  for leg_type in leg_types:
    try:
      leg_type_options.index(leg_type)
    except ValueError:
      raise ValueError('Invalid leg type.  Enter one or more of the following values for the leg type: {}'.format(', '.join(flight_options)))
      break
    leg_types[n]=''.join(['    ',leg_type])
    n=n+1

  n=0
  for flight in flights:
    try:
      flight_options.index(flight)
    except ValueError:
      raise ValueError('Invalid flight number.  Enter one of more of the following values for the flight number: {}'.format(', '.join(flight_options)))
    if flight=='all':
      flights=[' rf01',' rf02',' rf03',' rf04',' rf05',' rf06',' rf07',' rf08',' rf09',' rf10','rf11',' rf12',' rf13',' rf14',' rf15']
      break
    flights[n]=''.join([' ',flight])
    n=n+1
    
  #Load the flight table
  leg_table=np.loadtxt(path+'flight_table/flight_leg_table.txt',delimiter=',	',dtype=object)

  #Select legs that match all of the input requirements
  leg_length=pd.to_numeric(leg_table[:,5])-pd.to_numeric(leg_table[:,4])
  leg_type_index=np.where(np.isin(leg_table[:,3],leg_types) == True)
  time_index=np.where(leg_length >= min_time)
  flight_index=np.where(np.isin(leg_table[:,1],flights) == True)
  subset=np.intersect1d(np.intersect1d(leg_type_index,time_index),flight_index)

  #Print leg data to the terminal
  for subset_ind in subset:
    print('Flight: {}'.format(''.join(leg_table[subset_ind,1])))
    print('Leg Type: {}'.format(''.join(leg_table[subset_ind,3])))
    print('Index: {}'.format(''.join(leg_table[subset_ind,4])),' - {}'.format(''.join(leg_table[subset_ind,5])),)
    print('Time: {}'.format(''.join(leg_table[subset_ind,6])),' - {}'.format(''.join(leg_table[subset_ind,7])),)
    print('Latitude: {}'.format(''.join(leg_table[subset_ind,8])),' - {}'.format(''.join(leg_table[subset_ind,9])),)
    print('Longitude: {}'.format(''.join(leg_table[subset_ind,10])),' - {}'.format(''.join(leg_table[subset_ind,11])),)
    print('Altitude: {}'.format(''.join(leg_table[subset_ind,12])),' - {}'.format(''.join(leg_table[subset_ind,13])),)
    print('Temperature: {}'.format(''.join(leg_table[subset_ind,14])),' - {}'.format(''.join(leg_table[subset_ind,15])),)
    print('')
    print('')

  #If MakeFile=True, make a new text file with the subsetted leg data
  if MakeFile:

    header='Leg types: A = above cloud level leg, B = below cloud level leg, C = in-cloud level leg, U = Upward sounding, D = Downward sounding, S = Sawtooth, F = Ferry, T = Tarmac, CL = Clear.     Latitude, longitude and altitude are the GPS references (GGLAT,GGLON<GGALT) and temperature is the ambient reference (ATX).     Columns: leg number (whole campaign), flight number, leg number (flight), leg type, index (start), index (end), time (start), time (end), latitude (start), latitude (end), latitude(mean), longitude (start), longitude (end), longitude(mean), altitude (min), altitude (max), altitude(mean), temperature (min), temperature (max), temperature (mean)'

    np.savetxt(os.getcwd()+'/flight_leg_table_subset.txt',leg_table[subset,:],fmt='%5s',delimiter=',	',header=header)

  #If MakeArray=True, return an array with the subsetted data
  if MakeArray:
     print(leg_table[subset,:])
     return leg_table[subset,:]

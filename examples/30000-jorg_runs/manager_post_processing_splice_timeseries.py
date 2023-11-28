# Import required modules
import json
import os
import sys
import glob
import subprocess
import json

import numpy as np
from pymech.neksuite import field


list_with_case_names = sorted(glob.glob("*-rbc_*"))
case_name = list_with_case_names
number_of_cases = len(case_name)

for i in range(0,number_of_cases):
    case=(case_name[i]).strip()
    print("----------------------------------------------------")
    print("----------------------------------------------------")
    print("-------------Current case to check: " +case+" ------ ")
    print("----------------------------------------------------")

    cont = input("Continue? [yes/no]")
    if cont=="yes":

            list_with_data_folders = sorted(glob.glob(case+"/"+"*data_*"))
            data_folders = list_with_data_folders
            number_of_data_folders= len(data_folders)
            
            if number_of_data_folders==0 : 
                print("This case has no data folders, check if the data is in the main one or run the case")
                continue

            print("Data folders in this case:")
            for folder in data_folders:
                data_folder_string = folder.strip()
                data_folder = data_folder_string[len(case)+1:]
                print(data_folder)

            print("splicing together the nusselt text files found")
            path_to_spliced_file = "timeSeries"


            # Collect the data from the multiple data paths
            ts_all = []
            time_series_length = 0
            for folder in data_folders:
                data_folder_string = folder.strip()
                data_folder = data_folder_string[len(case)+1:]
                
                # Read the nusselt data
                nusselt_file_name = "nusselt.txt"
                input_nusselt_file = case + "/" + data_folder + "/" + nusselt_file_name
                print("Reading nusselt_from: "+input_nusselt_file)
                data = np.loadtxt(input_nusselt_file)

                # Read the polynomial order
                list_with_fields = sorted(glob.glob(case + "/" + data_folder + "/" + "*0.f0*"))
                input_field_file = list_with_fields[0].strip()
                header = field.read_header(input_field_file)
                lx = header.orders[0]

                # Read the rayleigh number of the case
                f = open (case + "/rayleigh.case", "r") 
                params_file = json.loads(f.read())
                f.close()
                Ra = params_file["case"]["fluid"]["Ra"]

                extra_params = 2
                extra_params += 3 # Nu_eps and tke
 
                # Have the nusselt data but add lx and Ra
                ts = np.zeros((data.shape[0], data.shape[1]+extra_params))
                ts[:,:4] = data[:,:]
                ts[:,4] = np.ones((data.shape[0]))*lx
                ts[:,5] = np.ones((data.shape[0]))*Ra

                # Check if nu_eps.txt exist
                nu_eps_file_name = "nu_eps.txt"
                input_nu_eps_file = case + "/" + data_folder + "/" + nu_eps_file_name
                if os.path.exists(input_nu_eps_file): 
                    print('nu_eps exist')
                    data_eps = np.loadtxt(input_nu_eps_file)
                    ts[:,6] = data_eps[:,1]
                    ts[:,7] = data_eps[:,2]
                else:
                    print('Nu_eps does not exist')
                
                # Check if tke.txt exist
                tke_file_name = "tke.txt"
                input_tke_file = case + "/" + data_folder + "/" + tke_file_name
                if os.path.exists(input_tke_file): 
                    print('tke exist')
                    data_tke = np.loadtxt(input_tke_file)
                    ts[:,8] = data_tke[:,1]
                else:
                    print('tke does not exist')

                time_series_length += data.shape[0]

                ts_all.append(ts)

            # Splice the data together
            ts_spliced = np.zeros((1, data.shape[1]+extra_params))
            print(len(ts_all))
            for i in range(0,len(ts_all)):
                ts_spliced = np.append(ts_spliced, ts_all[i],axis=0)

                # Check if the time is always increasing, this could be broken if the run stopped suddenly and the time series was not cleaned. 
                increasing = np.all(np.diff(ts_spliced[:,0]) > 0) 
                if increasing == False: 
                    print("the time series is not strictly increasing, overwriting instances in which the time is larger than a subsequent simulation") 
                    index_where_time_goes_back  = np.where(np.diff(ts_spliced[:,0]) < 0) 
                    print(index_where_time_goes_back[0][0]) 
                    index_that_should_be_deleted = np.where(ts_spliced[:index_where_time_goes_back[0][0]+1,0] >= ts_spliced[index_where_time_goes_back[0][0]+1,0])[0] 
                    print("The following entries (rows) of the time series will be deleted to favor the newer simulation") 
                    print(index_that_should_be_deleted)

                    print(ts_spliced.shape)
                    print(len(index_that_should_be_deleted))
                    ts_spliced = np.delete(ts_spliced, index_that_should_be_deleted, axis=0)
                    print(ts_spliced.shape)

                # Check again if delete worked
                increasing = np.all(np.diff(ts_spliced[:,0]) > 0) 
                print("is the series strictly increasing now?")
                print(increasing)


            print(ts_spliced[0,:])
            ts_spliced = ts_spliced[1:,:]

            print("Data was spliced, now interpolating")

            # Find the min dt and interpolate to that
            dt = 10000
            for i in range(0,ts_spliced.shape[0]-1):
                dt = min(dt, abs(ts_spliced[i+1,0]-ts_spliced[i,0]))
            print(dt)

            regular_t = np.arange(ts_spliced[0,0],ts_spliced[-1,0],dt)
            ts_interpolated = np.zeros((regular_t.shape[0], data.shape[1]+extra_params))
            ts_interpolated[:,0] = regular_t[:]
            for i in range(1, data.shape[1]+extra_params):
                ts_interpolated[:,i] = np.interp(regular_t[:], ts_spliced[:,0], ts_spliced[:,i])


            print("data was interpolated, now writting")

            # Save the array:
            np.savetxt(path_to_spliced_file + "/"+ case + "_original.txt", ts_spliced)
            np.savetxt(path_to_spliced_file + "/"+ case + "_interpolated.txt", ts_interpolated)

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
            new_data_folders = []
            for folder in data_folders:
                data_folder_string = folder.strip()
                data_folder = data_folder_string[len(case)+1:]
                print(data_folder)
                print(data_folder[6:])

                new_suffix = data_folder[6:].zfill(5)
                new_folder_string= data_folder_string[:len(case)+1] + data_folder[:6] + new_suffix

                print("old folder: "+data_folder_string)
                print("new folder: "+new_folder_string)

                command = "mv "+data_folder_string+ " " + new_folder_string  + ";" 
                print(command)
                process = subprocess.run(command, capture_output=True, shell=True)
                print(process.stdout.decode())
                print("====== moved =====")



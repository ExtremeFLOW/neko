# Import required modules
import json
import os
import sys
import glob
import subprocess




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


        path = "./"+case+"/"

        submit = input("Do you want to submit now? [yes/no]")
        if submit == "yes":
            command = "cd "+path+";"
            command = command + "sbatch run.sh " 
            print(command)
            process = subprocess.run(command, capture_output=True, shell=True)
            print(process.stdout.decode())
            print("====== Submitted =====")


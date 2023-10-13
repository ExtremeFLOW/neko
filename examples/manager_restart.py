# Import required modules
import json
import os
import sys
import glob
import subprocess




#list_with_case_names = glob.glob("*-rbc_*")

list_with_case_names = open("case_list.txt", "r")
case_name = list_with_case_names.readlines()
number_of_cases = len(case_name)
list_with_case_names.close()

for i in range(0,number_of_cases):
    case=(case_name[i]).strip()
    print("Current case to restart: " +case)
    cont = input("Continue? [yes/no]")
    if cont=="yes":
        path = "./"+case+"/"

        #===================================================
        joblimit_files = sorted(glob.glob(path+"joblimit*"))
        if joblimit_files == []:
            print("there is no job limit file")
            joblimit_index = "no_file"
        else:
            joblimit_index = joblimit_files[-1]
            print("The last job limit file encoutered was: " + joblimit_index)
        #===================================================
        fluid_files = sorted(glob.glob(path+"fluid0*"))
        if fluid_files == []:
            print("there is no fluid checkpoint file")
            fluid_index = "no_file"
        else:
            fluid_index = fluid_files[-1]
            print("The last fluid file encoutered was: " + fluid_index)
        #===================================================
        lastcheck_files = sorted(glob.glob(path+"rbc*"))
        if lastcheck_files == []:
            print("there is no rbc check file file")
            lastcheck_index = "no_file"
        else:
            lastcheck_index = lastcheck_files[-1]
            print("The last rbc checkpoint file encoutered was: " + lastcheck_index)

        option = int(input("Choose from which file to restart: [1]: joblimit.chkp, [2]: fluid.chkp, [3]: last rbc chkp "))

        file_names = [joblimit_index, fluid_index,lastcheck_index]


        #===================================================
        fld_files = sorted(glob.glob(path+"field0.f*"))
        if fld_files == []:
            print("No fld file encountered")
            fld_index=int(input("set the index that you want the restart file to have: "))            
        else:
            fld_index = int(fld_files[-1][-5:])
            print("using fld index: " + repr(fld_index))

        # Set the name of the checkpoint file
        chkp_file = "rbc_t" +repr(fld_index)+".chkp" 

        checkpoint_updated = False

        print("proceeding to copy it to "+path+chkp_file)
        command = ""
        command = command + "mv "+file_names[option-1]+" "+path+chkp_file 
        print(command)
        process = subprocess.run(command, capture_output=True, shell=True)
        print(process.stdout.decode())
        
        checkpoint_updated = True

        print("Reading the JSON file")
        f = open (path+"rayleigh.case", "r") 
        params_file = json.loads(f.read())
        f.close()
       
        params_file["case"]["end_time"] = 700
        params_file["case"]["job_timelimit"] = "23:45:00"
        params_file["case"]["restart_file"] = chkp_file

        params_file_str = json.dumps(params_file,indent = 4)
        print(params_file_str)

        with open(path+"rayleigh.case", "w") as outfile:
            outfile.write(params_file_str)

        #clean the directoory
        if fld_files != []:
            print("Cleaning directory")
            command = "cd "+path+";"
            command = command + "mkdir data ;" 
            command = command + "./move.sh ;" 
            command = command + "mv data data_t"+repr(fld_index) 
            print(command)
            process = subprocess.run(command, capture_output=True, shell=True)
            print(process.stdout.decode())

        submit = input("Do you want to submit now? [yes/no]")
        if submit == "yes":
            command = "cd "+path+";"
            command = command + "sbatch run.sh " 
            print(command)
            process = subprocess.run(command, capture_output=True, shell=True)
            print(process.stdout.decode())
            print("====== Submitted =====")

    else:
        print("Checking next case")

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
    print("Current case to check: " +case)


    path = "./"+case+"/"

    crash_files = sorted(glob.glob(path+"log.e*"))
    if crash_files == []:
        print("there is no log file")
        crash_index = "no_file"
    else:
        crash_index = crash_files[-1]
        print("The last error log file encoutered was: " + crash_index)

    logfile = open(crash_index, "r")
    log = logfile.readlines()
    number_of_lines = len(log)
    logfile.close()

    if number_of_lines > 3:
        print("WARNING - CRASH ENCOUNTERED")
        seelog=input("do you want to see the log? [yes/no]")
        if seelog == "yes":
            for line in log:
                print(line)

    else:
        print("No crash, you can restart")


     
    #cont = input("Continue? [yes/no]")

    #if cont=="yes":

    #else:
    #    print("Checking next case")

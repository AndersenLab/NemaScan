#!/usr/bin/env python

import os
#import sys  
# Read in the command line arguments

# The species of the strain set variants

#workflow_dir = sys.argv[1]
workflow_dir = "/projects/b1059/projects/Ryan/ortholog_sims/NemaScan"
#sp = sys.argv[2]
sp = "c_elegans"
#mapping_panel = sys.argv[3]
mapping_panel = "underground.gartersnake"

# Check if the workflow directory exists
dir_name = workflow_dir + "/test_data/" + sp + "/" + mapping_panel + "/plink_files"
#print("Checking " + dir_name)
cache_check = os.path.exists(dir_name)
#print(cache_check)

#Check if all the files are there if any of them are missing end the loop and call cache_check = False
plink_file_list = ["TO_SIMS.bed", "TO_SIMS.bim", "TO_SIMS.fam", "TO_SIMS.log", "TO_SIMS.nosex", "TO_SIMS.ped", "TO_SIMS.map"]
for file in plink_file_list:
    file_name = dir_name + "/" + file
    #print("Checking " + file_name)
    cache_check = os.path.exists(file_name)
    #print(cache_check)
    if cache_check == False:
        break

print(cache_check)

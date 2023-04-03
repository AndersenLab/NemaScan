#!/usr/bin/env python3
import csv
import time

# For a given number of genes in an orthogroup, classify the structure
def classify_structure(n_sp):
    if n_sp == 0:
        c_sp = "none"
    if n_sp == 1:
        c_sp = "one"
    if n_sp > 1:
        c_sp = "many"
    return(c_sp)

# For a given orthogroup file, classify the structure of each orthogroup 
# and return a dictionary of orthogroups by id and structure
def classify_orthogroups(sp1, sp2, sp3, og_file):
    orthogroups_dict = {}
    with open(og_file, 'r') as orthogroups:
        for line in orthogroups:
            og_id = line.rstrip("\n").split(": ")[0]
            seqs = line.rstrip("\n").split(": ")[1].split(" ")
            #Count the number of genes for each species
            n_sp1 = 0
            n_sp2 = 0
            n_sp3 = 0
            for seq in seqs: 
                if seq.startswith(sp1):
                    n_sp1 = n_sp1 + 1
                if seq.startswith(sp2):
                    n_sp2 = n_sp2 + 1
                if seq.startswith(sp3):
                    n_sp3 = n_sp3 + 1
            #Classify structure
            c_sp1 = classify_structure(n_sp1)
            c_sp2 = classify_structure(n_sp2)
            c_sp3 = classify_structure(n_sp3)
            og_struct = c_sp1 + ":" + c_sp2 + ":" + c_sp3
            orthogroups_dict[og_id] = {"id":og_id, sp1:n_sp1, sp2:n_sp2, sp3:n_sp3, "structure":og_struct}
    return(orthogroups_dict)

def filter_orthogroups(og_dict, structure):
    #Filters the orthogroups by structure (e.g. "one:none:none")
    filtered_ogs = {}
    for og in og_dict.values():
        if og["structure"] == structure:
            filtered_ogs[og["id"]] = og
    return(filtered_ogs)

if __name__ == "__main__":
        
    #Define the path to the orthogroups file
    orthogroups = "input_data/all_species/orthogroups/02.21.22_orthogroups/Orthogroups.txt"
    
    #Define the transcript prefixes for each species used in GFF 
    sp1 = "Transcript" #elegans
    sp2 = "QX1410"
    sp3 = "transcript_CTROP"

    #Classify the orthogroups
    og_dicts = classify_orthogroups(sp1, sp2, sp3, orthogroups)

    #Filter the orthogroups by structure
    structure = "one:one:one"

    filtered_ogs = filter_orthogroups(og_dicts, structure)

    #Convert the dictionary to a csv file
    #create a file name with today's date
    ## save todays date as a string
    today = time.strftime("%Y%m%d-%H%M%S")

    file_name = today + "_orthogroup_structures" + ".csv"
    col_name=["id", sp1, sp2, sp3, "structure"]
    with open(file_name, 'w') as csvFile: #Change the name of the output file here
            wr = csv.DictWriter(csvFile, fieldnames=col_name)
            wr.writeheader()
            for ele in og_dicts.values():
                wr.writerow(ele)
    
    #Save the filtered orthogroups to a file
    file_name = today + "_filtered_orthogroups" + ".csv"
    col_name=["id", sp1, sp2, sp3, "structure"]
    with open(file_name, 'w') as csvFile: #Change the name of the output file here
            wr = csv.DictWriter(csvFile, fieldnames=col_name)
            wr.writeheader()
            for ele in filtered_ogs.values():
                wr.writerow(ele)
    
    # Save the filtered orthogroups to a list of IDs
    file_name = today + "_filtered_orthogroup_ids" + ".txt"
    with open(file_name, 'w') as txtFile: #Change the name of the output file here
        for ele in filtered_ogs.values():
            txtFile.write(ele["id"] + "\n")
            
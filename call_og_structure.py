#!/usr/bin/env python3
import csv

def classify_structure(n_sp):
    if n_sp == 0:
        c_sp = "none"
    if n_sp == 1:
        c_sp = "one"
    if n_sp > 1:
        c_sp = "many"
    return(c_sp)

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


orthogroups = "Orthogroups.txt"
sp1 = "Transcript" #elegans
sp2 = "QX1410"
sp3 = "transcript_CTROP"

og_dicts = classify_orthogroups(sp1, sp2, sp3, orthogroups)


col_name=["id", sp1, sp2, sp3, "structure"]
with open("export.csv", 'w') as csvFile:
        wr = csv.DictWriter(csvFile, fieldnames=col_name)
        wr.writeheader()
        for ele in og_dicts.values():
            wr.writerow(ele)
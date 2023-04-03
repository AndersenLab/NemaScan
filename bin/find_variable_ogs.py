import pandas as pd 
import random

# Load the strain set variants 
def load_strain_set_variants(strain_set_variant_file):
    """
    Load the strain set variants from a plink .bim file
        ** THIS IS ALSO DEFINED IN sim_og_effects.py **

    """
    #Read in only the first 2 columns of the bim file with tab delimiter
    strain_set_variants = pd.read_csv(strain_set_variant_file, sep="\t", header=None, usecols=[0,1], names=["CHROM", "ID"])

    return strain_set_variants
def load_passing_ogs(passing_ogs_file):
    """
    Load the passing orthogroups from a file
    """
    #Read lines into a list
    with open(passing_ogs_file) as f:
        passing_ogs = f.readlines()
        #Remove whitespace characters like `\n` at the end of each line
        passing_ogs = [x.strip() for x in passing_ogs]
    return passing_ogs
    
    


    return passing_ogs

def add_var_orthogroups(strain_set_variants, sp, passing_ogs):
    """
    ** THIS IS ALSO DEFINED IN sim_og_effects.py **

    Add the orthogroup information to the strain set variants

    Parameters
    ----------
    strain_set_variants : pandas dataframe
        The strain set variants with the following columns:
        0: chromosome
        1: variant id
    sp : string
        The species of the strain set variants
        Reads in the master snps file with columns:
        0: variant id (eg. 2:1000000)
        1: transcript id (eg. C5641.2)
        2: orthogroup id (eg. OG0000000)
    """
# load the master snps file for the species - HARD CODED FOR NOW
    if sp == "c_elegans":
        master_snps_file = "input_data/c_elegans/orthogroup_snps/02.21.22_orthogroups/celeg-master.snps.genes"
        print(master_snps_file)
    elif sp == "c_briggsae":
        master_snps_file = "input_data/c_briggsae/orthogroup_snps/02.21.22_orthogroups/cbrig-master.snps.genes"
    elif sp == "c_tropicalis":
        master_snps_file = "input_data/c_tropicalis/orthogroup_snps/02.21.22_orthogroups/ctrop-master.snps.genes"
    else:
        print("Species not recognized")
        return
    #Read in the master snps file
    master_snps = pd.read_csv(master_snps_file, header=None, names=['ID', "GENE", "ORTHOGROUP"], delimiter= ' ')

    
    #Filter the master snps file to only include passing orthogroups
    master_snps = master_snps[master_snps["ORTHOGROUP"].isin(passing_ogs)]

    #Left join the master snps file to the strain set variants so that the orthogroup ids are added to the strain set variants using the variant id as the key
    strain_set_variants = pd.merge(strain_set_variants, master_snps, how="left", left_on="ID", right_on="ID")

    #remove the variants that do not have orthogroup ids
    strain_set_variants = strain_set_variants.dropna()
    #return the strain set variants with the orthogroup ids
    return strain_set_variants

def get_var_ogs(strain_set_variants):
    """
    Get orthogroups that have variants in the strain set

    Parameters
    ----------
    strain_set_variants : pandas dataframe
        The strain set variants with the following columns:
        0: chromosome
        1: variant id
        2: transcript id
        3: orthogroup id
    """
    #Get the orthogroup ids for each variant
    var_ogs = strain_set_variants["ORTHOGROUP"].unique()

    return var_ogs

def sample_casual_ogs(all_var_ogs, num_ogs, n_reps):
    """
    Sample the casual orthogroups

    Parameters
    ----------
    all_var_ogs : list
        The orthogroups that have variants in all three species
    num_ogs : int
        The number of orthogroups to sample
    """
    #Sample the number of causal orthogroups and save them as a list. Repeat this n_reps times and save the results in a list
    causal_ogs = [random.sample(all_var_ogs, num_ogs) for i in range(n_reps)]
    return causal_ogs

    

    return casual_ogs

def write_causal_ogs_file(causal_ogs, causal_ogs_file):
    """
    Write the causal orthogroups to a file where each line is a list of causal orthogroups, each
    separated by a comma eg. OG0000000, OG0000001, OG0000002 each list of orthogroups is denoted by a new line
    Parameters
    ----------
    causal_ogs : list
        A list containing lists of causal orthogroups
    causal_ogs_file : string
        The file to write the causal orthogroups to
    """
    #Open the file to write to
    with open(causal_ogs_file, 'w') as f:
        #Loop through each list of causal orthogroups
        for causal_og_list in causal_ogs:
            #Write the list of causal orthogroups to the file
            f.write(",".join(causal_og_list) + "\n")
    return


if __name__ == "__main__":
    n_causal_ogs = 5
    n_reps = 10

    ce_vars = load_strain_set_variants("test_data/plink_files/TO_SIMS.bim")
    cb_vars = load_strain_set_variants("test_data/plink_files/TO_SIMS.bim") #FOR TESTING THEY ARE ALL THE SAME ATM
    ct_vars = load_strain_set_variants("test_data/plink_files/TO_SIMS.bim")
   
    #Load orthogroups that passed the structure filtering
    passing_ogs = load_passing_ogs("20230331-130654_filtered_orthogroup_ids.txt")

    #Add the orthogroup ids to the strain set variants  
    ce_vars = add_var_orthogroups(ce_vars, "c_elegans", passing_ogs)
    cb_vars = add_var_orthogroups(cb_vars, "c_elegans", passing_ogs) #FOR TESTING THEY ARE ALL THE SAME ATM
    ct_vars = add_var_orthogroups(ct_vars, "c_elegans", passing_ogs)

    #Get the orthogroups that have variants in the strain set
    ce_var_ogs = get_var_ogs(ce_vars)
    cb_var_ogs = get_var_ogs(cb_vars)
    ct_var_ogs = get_var_ogs(ct_vars)

    #Get the orthogroups that have variants in all three species 
    all_var_ogs = set(ce_var_ogs).intersection(cb_var_ogs, ct_var_ogs)
    
    #Convert the set to a list
    all_var_ogs = list(all_var_ogs)

    
    #Sample the casual orthogroups
    casual_ogs = sample_casual_ogs(all_var_ogs, n_causal_ogs, n_reps)

    #Write out a file of selected orthogroups where each line is a set of orthogroups
    write_causal_ogs_file(casual_ogs, "causal_ogs.txt")
    
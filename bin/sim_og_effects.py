import pandas as pd
import numpy as np
import sys

# Load the strain set variants 
def load_strain_set_variants(strain_set_variant_file):
    """
    Load the strain set variants from a plink .bim file
    ** THIS IS ALSO DEFINED IN simulate_orthogroup_vars.py **
    """
    #Read in only the first 2 columns of the bim file with tab delimiter
    strain_set_variants = pd.read_csv(strain_set_variant_file, sep="\t", header=None, usecols=[0,1], names=["CHROM", "ID"])

    return strain_set_variants

def add_var_orthogroups(strain_set_variants, sp, passing_ogs, master_snps_dir):
    """
    ** THIS IS ALSO DEFINED IN simulate_orthogroup_vars.py **

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
# load the master snps file for the species - the snp files are pulled from the defined master_snps_dir in the NF pipeline
    if sp == "c_elegans":
        master_snps_file = master_snps_dir + "/celeg-master.snps.genes"
        print(master_snps_file)
    elif sp == "c_briggsae":
        master_snps_file = master_snps_dir + "/cbrig-master.snps.genes"
    elif sp == "c_tropicalis":
        master_snps_file = master_snps_dir + "/ctrop-master.snps.genes"
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


def select_og_variants(og_variants, n_var_og = 1):
    """
    Group variants by orthogroup and select n_var_og variants from each orthogroup
    """
    #Group variants by orthogroup
    og_variants = og_variants.groupby("ORTHOGROUP")
    #Select n_var_og variants from each orthogroup
    og_variants = og_variants.apply(lambda x: x.sample(n_var_og))
    #Retrun a dataframe of the selected variants
    return og_variants

def simulate_og_effect_gamma(og_variants, n_var, og_effect_shape = 0.4, og_effect_scale = 1.66):
    effects = np.random.default_rng().gamma(og_effect_shape, og_effect_scale, n_var)
    directions = np.random.choice([-1,1], n_var)
    effects = effects * directions
    og_variants["EFFECT"] = effects
    return(og_variants)


if __name__ == "__main__":

    #Define orthogroups from command line arguments
    og1 = sys.argv[1]
    og2 = sys.argv[2]
    og3 = sys.argv[3]
    og4 = sys.argv[4]
    og5 = sys.argv[5]

    #Get the list of variants in the strain sets
    strain_set_variant_file = sys.argv[6]

    #Get the directory of the master snps files
    master_snps_dir = sys.argv[7]

    #Define the species
    sp = sys.argv[8]

    #og1 = "OG0010644"
    #og2 = "OG0010836"
    #og3 = "OG0011342"
    #og4 = "OG0011339"
    #og5 = "OG0010147"


    sim_ogs = [og1, og2, og3, og4, og5]

    #Read in annotated strain_set variants
    strain_var = load_strain_set_variants(strain_set_variant_file)


    #Add the passing orthogroup ids to the strain set variants and only return the variants that have orthogroup ids from the joining
    print("Adding orthogroup ids to strain set variants")
    og_vars = add_var_orthogroups(strain_var, sp, sim_ogs, master_snps_dir)
    
    #Select Causal variants for orthogroups
    causal_og_vars = select_og_variants(og_vars, 1)
    print("Selected Causal Variants")
    print(causal_og_vars)

    #Simulate effects for causal variants
    causal_og_vars = simulate_og_effect_gamma(causal_og_vars, 5)
    
    #Write output for trait simulations - just id and effect
    causal_og_vars[["ID", "EFFECT"]].to_csv("causal_og_vars.txt", sep = " ", index = False, header = False)
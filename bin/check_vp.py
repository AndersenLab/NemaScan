import argparse

#Parse command line arguments
def parse_commandline(): 
    """Parse the arguments given on the command-line.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check_vp",
                        help = "Path to the check_vp file")
    parser.add_argument("--simulated_phenos",
                        help = "path to the simulated phenotypes file")                 
    args = parser.parse_args()

    return args

def get_vp(check_vp_file):
    """
    # Read in the check_vp file a tab delim file with the following columns:
    # 1. Source
    # 2. Variance
    # 3. Standard error
    # And pull out variance for the Vp source
    """
    with open(check_vp_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0] == 'Vp':
                vp = line[1]
    return vp
def increase_pheno_var(phenos_file):
    """
    Read in the simulated phenotypes file a space delim file with the following columns:
    - Strain 1 
    - Strain 1 
    - Phenotype value
    Multiply the phenotype value by 1000 and save `new_phenos.temp`
    """
    with open(phenos_file, 'r') as f:
        for line in f:
            line = line.strip().split(' ')
            pheno = float(line[2])
            new_pheno = pheno * 1000
           # print(new_pheno)
            with open('new_phenos.temp', 'a') as f:
                f.write(line[0] + ' ' + line[1] + ' ' + str(new_pheno) + '\n')



if __name__ == "__main__":
    #Parse command line arguments
    # Parse command line arguments
    args = parse_commandline()
    check_vp_file = args.check_vp
    phenos_file = args.simulated_phenos

    # Get the Vp from the check_vp file
    vp = get_vp(check_vp_file)
    print(vp)
    
    #If VP is less than 0.000001 then print message 
    if float(vp) <= 0.000001:
        print("Vp is less than 0.000001")
        #update the phenos file by multiplying the 3rd colum by 1000
        increase_pheno_var(phenos_file)
    else:
        print("Vp is greater than 0.000001")

        


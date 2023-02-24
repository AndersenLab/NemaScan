#! usr/bin/env nextflow


include{prepare_simulation_files} from './modules/simulations.nf'

//ce_vcf = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.bcsq.vcf.gz")
//ce_vcf_index = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.bcsq.vcf.gz.tbi")

//cb_vcf = Channel.fromPath("/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.bcsq.vcf.gz")
//cb_vcf_index = Channel.fromPath("/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.bcsq.vcf.gz.tbi")

//ct_vcf = Channel.fromPath("/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.bcsq.vcf.gz")
//ct_vcf_index = Channel.fromPath("/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.bcsq.vcf.gz.tbi")

rename_key = Channel.fromPath("/projects/b1059/projects/Ryan/telomere_variation/NemaScan/input_data/all_species/rename_chromosomes")
maf_file = Channel.fromPath("/projects/b1059/projects/Ryan/telomere_variation/NemaScan/input_data/all_species/simulate_maf.csv").splitCsv()


workflow{

File pop_file = new File("/projects/b1059/projects/Ryan/telomere_variation/NemaScan/test_data/test_orthogroup_samples.txt") ;

sp_ids = [["c_elegans", "/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.bcsq.vcf.gz", "/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.bcsq.vcf.gz.tbi"], ["c_briggsae", "/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.bcsq.vcf.gz", "/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.bcsq.vcf.gz.tbi"], ["c_tropicalis", "/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.bcsq.vcf.gz", "/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.bcsq.vcf.gz.tbi"]]

Channel.from(pop_file.collect { it.tokenize( ' ' ) })
          .map {SP, SM, STRAINS -> [SP, SM, STRAINS] }
          .join(Channel.from(sp_ids), by:[0])
          .combine(rename_key)
          .combine(maf_file)
          .map {
            tuple ->[tuple[0], //extract species ID
                    (tuple[1]), //extract sample name
                    (tuple[2]), //extract sample list
                    file(tuple[3]), // convert path to file obj for vcf
                    file(tuple[4]), // index
                    file(tuple[5]), // rename key
                    tuple[6]] // MAF
        } | prepare_simulation_files 
        

}
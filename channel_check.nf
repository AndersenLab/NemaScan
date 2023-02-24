#!/ usr/bin/env nextflow


File pop_file = new File("/projects/b1059/projects/Ryan/telomere_variation/NemaScan/test_data/test_orthogroup_samples.txt")

// ogs = file("/projects/b1059/projects/Ryan/telomere_variation/NemaScan/selected_ogs.txt").readLines()

sp_ids = [["c_elegans", "path_1"] ,["c_briggsae", "path_2"], ["c_tropicalis", "path"]]


strain_its= Channel.from(pop_file.collect { it.tokenize( ' ' ) })
          .map {SP, SM, STRAINS -> [SP, SM, STRAINS] }
          .join(Channel.from(sp_ids), by:[0])

strain_its.view()






///id_sp_channel = Channel.from(ogs).combine(Channel.from(sp_ids))
//id_sp_channel.view()


// both = ogs.map{it -> [OG].combine(sp)}

//both.view()


process create_fasta_promoter {
    cpus = 8
    memory = '64 GB'
    time = '4h'

    label 'bedtools'

    input:
    tuple val(species), val(mt_contig),  file("${species}.gtf"), file("${species}.fa")

    output:
    tuple val(species), file("${species}_5000up_1000down_TSS.fasta")    

    script:
    """
    #!/bin/bash

  
    bedtools getfasta \
	-fi ${species}.fa \
	-bed ${projectDir}/data/TSSs/${species}_TSS_regions.bed \
	-fo ${species}_5000up_1000down_TSS.fasta -nameOnly


    """
}

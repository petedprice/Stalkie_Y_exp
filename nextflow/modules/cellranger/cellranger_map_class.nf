process cellranger_map_class {
    cache 'lenient' 
    cpus = 16
    memory = '64 GB'
    time = '4h'

    label 'samtoolsetc'

    tag {'Cellranger_count_' + '_' + species + '_' + sample }

    publishDir 'read_type_mapping', mode: 'copy', overwrite: true, pattern: '*read_types*'

    input:
    tuple val(species), val(sample), file("${sample}_crdata")

    output:
    tuple val(species), val(sample), file("${sample}_read_types.tsv")

    script:
    """
    #!/bin/bash
    python ${projectDir}/pyscripts/get_read_class.py \
	--bam ${sample}_crdata/outs/possorted_genome_bam.bam \
	--out ${sample}_read_types.tsv

    """

}

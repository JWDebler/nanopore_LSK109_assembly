
def helpMessage() {
    log.info"""
    # Nanopore LSK109 chemistry r9 flowcell genome assembly and polishing
    A pipeline for assembly and polishing of fungal genomes from Oxford Nanopore reads

    ## Examples
    nextflow run jwdebler/nanopore_LSK109_assembly -resume -latest -profile docker,nimbus --reads "reads/"

    ## Parameters
    --reads <glob>
        Required
        A folder containing 1 files per sample. 
        The basename of the file is used as the sample ID.
       
        Example of file names: `Sample1.fastq.gz`, `Sample2.fastq.gz`.
        (Default: a folder called `reads/`)

    --genomeSize <glob>
        not required
        Size of genome, for example "42m" (Default: 42m)

    --medakaModel <glob>
        not required
        Which basecaller model was used?
        r941_min_sup_g507 (kit109, sup)
        (Default: r941_min_sup_g507)

    --minlen
        Min read length to keep for assembly
        (Default: 1000)

    --quality
        Min read q-score to keep for read filtering
        (Default: 10)

    --canuSlow
        Disables canu fast mode.
        (Default: false)

    --outdir <path>
        The directory to store the results in.
        (Default: `assembly`)

    ## Exit codes
    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """
}

//by default script looks for reads in a directory called "reads"

params.reads="reads/"
params.size="42m"
params.medakaModel="r941_min_sup_g507"
params.quality="10"
params.minlen="1000"

if (params.help) {
    helpMessage()
    exit 0
}

nanoporeReads = Channel
    .fromPath(params.reads + "*.fast{q,q.gz}", checkIfExists: true)
    .map {file -> [file.simpleName, file]}
    .tap { ReadsForQC }
    .tap { ReadsForNanoplot }
    .view()

process version_canu {

    label "canu"

    output:
    path 'versions.txt' into canu_version

    """
    echo canu: >> versions.txt
    canu --version >> versions.txt
    echo --------------- >> versions.txt
    """
}


process version_medaka {

    label "medaka"

    output:
    path 'versions.txt' into medaka_version

    """
    echo medaka: >> versions.txt
    medaka --version >> versions.txt
    echo --------------- >> versions.txt
    """
}


process version_seqkit {

    label "seqkit"

    output:
    path 'versions.txt' into seqkit_version

    """
    echo seqkit: >> versions.txt
    seqkit version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process version_flye{

    label "flye"

    output:
    path 'versions.txt' into flye_version

    """
    echo flye: >> versions.txt
    flye --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process version_nextdenovo {

    label "nextdenovo"

    output:
    path 'versions.txt' into nextdenovo_version

    """
    echo nextdenovo: >> versions.txt
    nextDenovo --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process version_chopper {

    label "chopper"

    output:
    path 'versions.txt' into chopper_version

    """
    echo chopper: >> versions.txt
    chopper --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process version_minimap2 {

    label "chopper"

    output:
    path 'versions.txt' into minimap2_version

    """
    echo minimap2: >> versions.txt
    minimap2 --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process version_samtools {

    label "chopper"

    output:
    path 'versions.txt' into samtools_version

    """
    echo samtools: >> versions.txt
    samtools version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process versions {

    input:
    path "canu.txt" from canu_version
    path "medaka.txt" from medaka_version
    path "seqkit.txt" from seqkit_version
    path "flye.txt" from flye_version
    path "nextdenovo.txt" from nextdenovo_version
    path "chopper.txt" from chopper_version
    path "minimap2.txt" from minimap2_version
    path "samtools.txt" from samtools_version

    publishDir "${params.outdir}/", mode: 'copy', pattern: 'versions.txt'

    output:
    path "versions.txt"

    script:
    """
    cat canu.txt medaka.txt seqkit.txt flye.txt nextdenovo.txt chopper.txt minimap2.txt samtools.txt > versions.txt
    """
}

process QC_read_removal_and_filtering {

    label "chopper"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.fastq.gz'

    input:
    tuple sampleID, "${sampleID}.fastq.gz" from ReadsForQC

    output:
    tuple sampleID, "${sampleID}.chopper.200bp.q${params.quality}.fastq.gz" into FilteredDuplex200
    tuple sampleID, "${sampleID}.chopper.${params.minlen}bp.q${params.quality}.fastq.gz" into FilteredDuplex1000
    

    """
    wget https://raw.githubusercontent.com/JWDebler/nanopore_kit14_assembly/main/data/DCS.fasta
    minimap2 -d dcs.mmi DCS.fasta
    seqkit rmdup -n ${sampleID}.fastq.gz | minimap2 -t "${task.cpus}" -ax map-ont dcs.mmi - | samtools view -O fastq -@ "${task.cpus}" - | chopper -t ${task.cpus}  -q ${params.quality} -l 200 | pigz -9 > ${sampleID}.chopper.200bp.q${params.quality}.fastq.gz 
    chopper -i ${sampleID}.chopper.200bp.q${params.quality}.fastq.gz -t ${task.cpus} -l ${params.minlen} | pigz -9 > ${sampleID}.chopper.${params.minlen}bp.q${params.quality}.fastq.gz
    """
}



process QC_nanoplot_Raw {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, "${sampleID}.fastq.gz" from ReadsForNanoplot

    output:
    path "*.html"
    
    
    """
    NanoPlot \
    --fastq ${sampleID}.fastq.gz \
    -o output && \
    cp output/NanoPlot-report.html ${sampleID}.nanoplot.html
    """
}

process QC_nanoplot_Chopper {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, "${sampleID}.fastq.gz" from FilteredDuplex1000

    output:
    path "*.html"
    tuple sampleID, "${sampleID}.fastq.gz" into FilterdForAssembly
    
    """
    NanoPlot \
    --fastq ${sampleID}.fastq.gz \
    -o output && \
    cp output/NanoPlot-report.html ${sampleID}.nanoplot.chopper.html
    """
}

FilterdForAssembly
.tap { FilteredForFlye }
.tap { FilteredForNextdenovo }
.tap { FilteredForMedaka }

FilteredDuplex200
.tap { Reads200bpForMedakaFlye }
.tap { Reads200bpForMedakaNextdenovo }

// flye assembly
process Assembly_flye {

    label "flye"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/03-assembly", pattern: '*_flye.fasta'

    input:
    tuple sampleID, "reads.fastq.gz" from FilteredForFlye

    output:
    tuple sampleID, "${sampleID}_flye.fasta" into MedakaFlye

    """
     flye \
    --nano-hq reads.fastq.gz \
    --genome-size ${params.size} \
    --asm-coverage 50 \
    --threads "${task.cpus}" \
    --out-dir ${sampleID}.flye 

    cp ${sampleID}.flye/assembly.fasta ${sampleID}_flye.fasta
    """
}

process Assembly_nextdenovo {

    label "nextdenovo"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/03-assembly", pattern: '*_nextdenovo.fasta'
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*corredted.fasta'

    input:
    tuple sampleID, "reads.fastq.gz" from FilteredForNextdenovo

    output:
    tuple sampleID, "${sampleID}_nextdenovo.fasta" into MedakaNextdenovo
    tuple sampleID, "${sampleID}.nextdenovo.corredted.fasta"

    """
    ls reads.fastq.gz > ${sampleID}.fofn

    echo '''
    [General]
    job_type = local
    job_prefix = ${sampleID}.nextdenovo
    task = all
    rewrite = yes
    deltmp = yes
    parallel_jobs = 5
    input_type = raw
    read_type = ont # clr, ont, hifi
    input_fofn = ${sampleID}.fofn
    workdir = ${sampleID}.nextdenovo

    [correct_option]
    read_cutoff = 1k
    genome_size = ${params.size} # estimated genome size
    sort_options = -m 4g -t 3
    minimap2_options_raw = -t 3
    pa_correction = 5
    correction_options = -p 10

    [assemble_option]
    minimap2_options_cns = -t 8
    nextgraph_options = -a 1
    ''' > ${sampleID}.config

    nextDenovo ${sampleID}.config

    cp ${sampleID}.nextdenovo/03.ctg_graph/nd.asm.fasta ${sampleID}_nextdenovo.fasta

    cat ${sampleID}.nextdenovo/02.cns_align/01.seed_cns.sh.work/seed_cns*/cns.fasta > ${sampleID}.nextdenovo.corredted.fasta
    """
}

process Polishing_medaka_flye {

    label "medaka"
    tag {sampleID}

    input:
    tuple sampleID, "${sampleID}.flye.fasta", "${sampleID}.fastq.gz" from MedakaFlye.join(Reads200bpForMedakaFlye)

    output:
    tuple sampleID, "${sampleID}_flye_medaka.fasta" into SeqkitFlye

    """
    
    medaka_consensus \
    -d flye.fasta \
    -i reads.fastq.gz \
    -o ${sampleID}_medaka_output \
    -t ${task.cpus} \
    -m ${params.medakaModel}

    cp ${sampleID}_medaka_output/consensus.fasta ${sampleID}_flye_medaka.fasta
    """
}

process Polishing_medaka_nextdenovo {

    label "medaka"
    tag {sampleID}

    input:
    tuple sampleID, "${sampleID}.nextdenovo.fasta", "${sampleID}.fastq.gz"  from MedakaNextdenovo.join(Reads200bpForMedakaNextdenovo)

    output:
    tuple sampleID, "${sampleID}_nextenovo_medaka.fasta" into SeqkitNextdenovo

    """
    medaka_consensus \
    -d nextdenovo.fasta \
    -i reads.fastq.gz \
    -o ${sampleID}_medaka_output \
    -t ${task.cpus} \
    -m ${params.medakaModel}

    cp ${sampleID}_medaka_output/consensus.fasta ${sampleID}_nextenovo_medaka.fasta
    """
}

process Cleanup_seqkitFlye {

    label "seqkit"
    tag {sampleID}

    publishDir "${params.outdir}/${sampleID}/04-medaka-polished", pattern: '*.fasta'

    input:
    tuple sampleID, "${sampleID}_flye_medaka.unsorted.fasta" from SeqkitFlye

    output:
    tuple sampleID, "${sampleID}_flye_medaka.fasta" into FlyeForRagtag

    """
    seqkit sort -lr ${sampleID}_flye_medaka.unsorted.fasta > ${sampleID}_flye_medaka.sorted.fasta
    seqkit replace -p '.+' -r 'flye_ctg_{nr}' --nr-width 2 ${sampleID}_flye_medaka.sorted.fasta > ${sampleID}_flye_medaka.fasta
    """
}

process Cleanup_seqkitNextdenovo {

    label "seqkit"
    tag {sampleID}

    publishDir "${params.outdir}/${sampleID}/04-medaka-polished", pattern: '*.fasta'

    input:
    tuple sampleID, "${sampleID}_nextdenovo_medaka.unsorted.fasta" from SeqkitNextdenovo

    output:
    tuple sampleID, "${sampleID}_nextdenovo_medaka.fasta" into NextDenovoForRagtag


    """
    seqkit sort -lr ${sampleID}_nextdenovo_medaka.unsorted.fasta > ${sampleID}_nextdenovo_medaka.sorted.fasta
    seqkit replace -p '.+' -r 'nd_ctg_{nr}' --nr-width 2 ${sampleID}_nextdenovo_medaka.sorted.fasta > ${sampleID}_nextdenovo_medaka.fasta
    """
}

// compare assemblies with RagTag and order according to Nextdenovo

process Cleanup_ragtag {

    label "ragtag"
    tag {sampleID}

    publishDir "${params.outdir}/${sampleID}/05-ragtag"

    input:
    tuple sampleID, "nextdenovo.fasta", "flye.fasta" from NextDenovoForRagtag.join(FlyeForRagtag)

    output:
    path "ragtag.*"

    """
    ragtag.py scaffold nextdenovo.fasta flye.fasta
    cp ragtag_output/* .
    """
}

// // read correction
// process Correction_canu {

//     label "canu"
//     tag {sampleID}
//     publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.fasta.gz'
//     publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.report'

//     input:
//     tuple sampleID, "reads.fastq" from ReadsForCorrection

//     output:
//     path "${sampleID}.corrected.fasta.gz"
//     path "${sampleID}.corrected.report"

//     script:
//     // See: https://groovy-lang.org/operators.html#_elvis_operator
//     fast_option = params.canuSlow ? "" : "-fast "

//     """
//     canu \
//     -correct \
//     -p ${sampleID} \
//     -d ${sampleID} \
//     genomeSize=${params.size} \
//     ${fast_option} \
//     -nanopore reads.fastq

//     cp ${sampleID}/*correctedReads.fasta.gz ${sampleID}.corrected.fasta.gz
//     cp ${sampleID}/*.report ${sampleID}.corrected.report
//     """
// }

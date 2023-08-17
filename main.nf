
// Files need to be duplex called and split into "sampleID.duplex.fastq.gz" and "sampleID.simplex.fastq.gz"


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

if (params.help) {
    helpMessage()
    exit 0
}

nanoporeReads = Channel
    .fromPath(params.reads, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap { ReadsForDCSQC }
    .view()

process canu_version {

    label "canu"

    output:
    path 'versions.txt' into canu_version

    """
    echo canu: >> versions.txt
    canu --version >> versions.txt
    echo --------------- >> versions.txt
    """
}


process medaka_version {

    label "medaka"

    output:
    path 'versions.txt' into medaka_version

    """
    echo medaka: >> versions.txt
    medaka --version >> versions.txt
    echo --------------- >> versions.txt
    """
}


process seqkit_version {

    label "seqkit"

    output:
    path 'versions.txt' into seqkit_version

    """
    echo seqkit: >> versions.txt
    seqkit version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process flye_version {

    label "flye"

    output:
    path 'versions.txt' into flye_version

    """
    echo flye: >> versions.txt
    flye --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process nextdenovo_version {

    label "nextdenovo"

    output:
    path 'versions.txt' into nextdenovo_version

    """
    echo nextdenovo: >> versions.txt
    nextDenovo --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process chopper_version {

    label "chopper"

    output:
    path 'versions.txt' into chopper_version

    """
    echo chopper: >> versions.txt
    chopper --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process minimap2_version {

    label "minimap2"

    output:
    path 'versions.txt' into minimap2_version

    """
    echo minimap2: >> versions.txt
    minimap2 --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process samtools_version {

    label "samtools"

    output:
    path 'versions.txt' into samtools_version

    """
    echo samtools: >> versions.txt
    samtools version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process version {

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

process minimap_DCS {

    label "minimap2"
    tag {sampleID}

    input:
    tuple sampleID, 'reads.fastq.gz' from ReadsForDCSQC

    output:
    tuple sampleID, 'reads.sam'into DCSalignments

    """
    wget https://raw.githubusercontent.com/JWDebler/nanopore_kit14_assembly/main/data/DCS.fasta
    minimap2 -d dcs.mmi DCS.fasta
    minimap2 -t "${task.cpus}" -ax map-ont dcs.mmi reads.fastq.gz > reads.sam
    """

}

process filter_DCS_reads {

    label "samtools"
    tag {sampleID}

    input:
    tuple sampleID, 'reads.sam' from DCSalignments

    output:
    tuple sampleID, 'reads.fastq' into DCSFilteredReads

    """
    samtools view -@ "${task.cpus}" -b -f 4 reads.sam | samtools fastq -@ "${task.cpus}" - > reads.fastq
    """
}

DCSFilteredReads
.tap { ReadsForChopper }
.tap { ReadsForQC }
.tap { ReadsForCorrection }

// filtering reads

process chopper {

    label "chopper"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.fastq.gz'

    input:
    tuple sampleID, 'reads.fastq' from ReadsForChopper

    output:
    path "${sampleID}.chopper.fastq.gz"
    tuple sampleID, "${sampleID}.chopper.fastq.gz" into FilteredReads

    """
    cat reads.fastq | chopper -q 10 -l 200 | gzip -9 > ${sampleID}.chopper.fastq.gz

    """
}


process nanoplot_Raw_Simplex {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, 'reads.fastq' from ReadsForQC

    output:
    path "*.html"
    
    
    """
    NanoPlot \
    --fastq reads.fastq \
    -o output && \
    cp output/NanoPlot-report.html ${sampleID}.nanoplot.html
    """
}

process nanoplot_Chopper_Simplex {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, 'reads.fastq.gz' from FilteredReads

    output:
    path "*.html"
    tuple sampleID, "reads.fastq.gz" into FilterdForAssembly
    
    """
    NanoPlot \
    --fastq reads.fastq.gz \
    -o output && \
    cp output/NanoPlot-report.html ${sampleID}.nanoplot.chopper.html
    """
}

FilterdForAssembly
.tap { FilteredForFlye }
.tap { FilteredForNextdenovo }
.tap { FilteredForMedaka }


// flye assembly
process flye {

    label "flye"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/03-assembly", pattern: '*_flye.fasta'

    input:
    tuple sampleID, "reads.fastq.gz", from FilteredForFlye

    output:
    tuple sampleID, "${sampleID}_flye.fasta", "reads.fastq.gz" into MedakaFlye

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

process nextdenovo {

    label "nextdenovo"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/03-assembly", pattern: '*_nextdenovo.fasta'
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*corredted.fasta'

    input:
    tuple sampleID, "reads.fastq.gz" from FilteredForNextdenovo

    output:
    tuple sampleID, "${sampleID}_nextdenovo.fasta", "reads.fastq.gz" into MedakaNextdenovo
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

process medaka_flye {

    label "medaka"
    tag {sampleID}

    input:
    tuple sampleID, "flye.fasta", "reads.fastq.gz" from MedakaFlye

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

process medaka_nextdenovo {

    label "medaka"
    tag {sampleID}

    input:
    tuple sampleID, "nextdenovo.fasta", "reads.fastq.gz"  from MedakaNextdenovo

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

process seqkitFlye {

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

process seqkitNextdenovo {

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

process ragtag {

    label "ragtag"
    tag {sampleID}

    publishDir "${params.outdir}/${sampleID}/", saveAs: '05-ragtag'

    input:
    tuple sampleID, "nextdenovo.fasta", "flye.fasta" from NextDenovoForRagtag.join(FlyeForRagtag)

    output:
    dir "ragtag_output"

    """
    ragtag.py scaffold nextdenovo.fasta flye.fasta
    """
}

// read correction
process canu {

    label "canu"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.fasta.gz'
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.report'

    input:
    tuple sampleID, "reads.fastq" from ReadsForCorrection

    output:
    path "${sampleID}.corrected.fasta.gz"
    path "${sampleID}.corrected.report"

    script:
    // See: https://groovy-lang.org/operators.html#_elvis_operator
    fast_option = params.canuSlow ? "" : "-fast "

    """
    canu \
    -correct \
    -p ${sampleID} \
    -d ${sampleID} \
    genomeSize=${params.size} \
    ${fast_option} \
    -nanopore reads.fastq

    cp ${sampleID}/*correctedReads.fasta.gz ${sampleID}.corrected.fasta.gz
    cp ${sampleID}/*.report ${sampleID}.corrected.report
    """
}

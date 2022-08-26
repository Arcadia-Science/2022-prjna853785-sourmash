#!/usr/bin/env nextflow 

nextflow.enable.dsl = 2

/* to include processes that are specified in a separate file/folder:
 include {p2p; send; sketch; taxonomy; receive_filter; receive_nofilter} from './modules/processes.nf'*/

params.fasta = "inputs/raw/*.fasta"
params.ksize  = 31 // set the default kmer size to 31. This can be controlled on the command line with nextflow run main.nf --ksize 21
// note I'm still calculating the sketches with all of the information so that this is always saved, but it will control the behavior of compare
// and later, it will also control the behavior of gather.
// also note these can be control with a parameters file like params.json https://carpentries-incubator.github.io/workflows-nextflow/03-workflow_parameters/index.html#parameter-file
// The parameter is referenced within a process shell script with ${ksize}


process sourmash_sketch {
    tag { sample } // define a tag that should carry sample information and control output file names...NOT SURE HOW THIS WORKS YET
    publishDir 'nextflow/sourmash_sketch' // files are published to the target folder creating a symbolic link for each process output that links the file produced into the process working directory. 
    // ${projectDir} is an implicit variable that specifies the directory where the main script is located.
    // We use it here bc nextflow scripts are executed in a separate working dir

    input: // defines the input dependencies, which determines the number of times a process is executed
    tuple val(sample_id), path(fasta)
    
    output: // defines the output channels used by the process
    path "${sample_id}.sig"

    script: // defines commands executed by the process
    /* By default, commands are run in bash.
       To change to a different language, a Shebang declaration is used.
       For example,  #!/usr/bin/env Rscript will turn the quote block into an R script.
       When the script gets to be long, it's recommended that the script is abstracted to its own file, and placed in the bin folder.
       This folder gets placed in the PATH so the script can be run. */
    """ 
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name "${sample_id}" -o "${sample_id}.sig" ${fasta}
    """
    // in the script block, bash variables need to be espaced with \ 
}

/*process sourmash_compare {
    input:
    
    output:
    
    shell:
    
}*/

// one workflow will be sourmash compare
workflow {
    // Make a channel containing fasta files. 
    // The channel will emit a tuple containing the file path and the file basename. 
    // The file basename will be used as a sample ID for output files and sig names.
    fasta_ch = Channel.fromPath(params.fasta).map {file -> tuple(file.baseName, file)}
    fasta_ch.view()
    sourmash_sketch(fasta_ch) // | sourmash_compare | sourmash_plot 
}


// you can have two channels as output. And then they will be 0 indexed.

/* TODO ~Monday:
+ figure out many:1 for sourmash_compare
+ add sourmash plot process
+ figure out docker containers for running, and document the build process
TODO longer term:
+ figure out how to re-use the sourmash sketch process -- or should there be two workflows, one that does sketch -> compare -> plot, and one that does sketch -> gather etc. Probably not two workflows bc I only want to have to sketch files once.
+ figure out how to encode download for dbs. Probably should parameterize and if/else on kmer size.


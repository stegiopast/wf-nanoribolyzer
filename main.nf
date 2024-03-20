#!/usr/bin/ nextflow
@Grab('com.xlson.groovycsv:groovycsv:1.1')
import static com.xlson.groovycsv.CsvParser.parseCsv
import java.io.File;
nextflow.enable.dsl=2

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                     Checkup imported Variables                                                          //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

println "Seqfolder"
println params.sample_folder
println ""

println "Color"
println params.color
println ""

println "Script folder"
println params.script_folder
println ""

println "Output folder"
println params.output_folder
println ""

println "Basecalling model"
println params.basecalling_model
println ""

println "Threads"
println params.threads
println ""

println "Sample type"
println params.sample_type
println ""


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                     Basecalling Dorado                                                                  //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process dorado_basecalling{
    label 'dorado_basecaller'
    publishDir "${params.out_dir}", mode: 'copy'
    input:
    path(sample_folder)
    path(basecalling_model)

    output:
    val done
    path('basecalling_output/basecalled.bam'), emit: basecalled_bam
    path('basecalling_output/sequencing_summary.txt')
    path('basecalling_output/basecalled_not_trimmed.fastq.gz'), emit: fastq_not_trimmed
    path('converted_to_pod5/converted.pod5'), emit: converted_pod5
    
    script:
    done = 1
    """ 
    mkdir -p basecalling_output
    mkdir -p converted_to_pod5
    (ls ${sample_folder}/*.pod5 && export filetype=pod5) || export filetype=fast5 
    if [ \$filetype == fast5 ]
    then 
        pod5 convert fast5 ${sample_folder}/*.fast5\
         --output converted_to_pod5/converted.pod5\
         --force-overwrite
        dorado basecaller ${basecalling_model} converted_to_pod5/\
         > basecalling_output/basecalled.bam 
        dorado summary basecalling_output/basecalled.bam\
         > basecalling_output/sequencing_summary.txt
        samtools bam2fq basecalling_output/basecalled.bam\
         -@ ${params.threads}\
         > basecalling_output/basecalled_not_trimmed.fastq
        gzip basecalling_output/basecalled_not_trimmed.fastq\
         -c\
         -1\
          > basecalling_output/basecalled_not_trimmed.fastq.gz
    fi
    if [ \$filetype == pod5 ]
    then
        dorado basecaller ${basecalling_model} ${sample_folder}/*.pod5\
         > basecalling_output/basecalled.bam
        dorado summary basecalling_output/basecalled.bam\
         > basecalling_output/sequencing_summary.txt
        samtools bam2fq basecalling_output/basecalled.bam\
         -@ ${params.threads}\
         > basecalling_output/basecalled_not_trimmed.fastq
        gzip basecalling_output/basecalled_not_trimmed.fastq\
         -c\
         -1\
          > basecalling_output/basecalled_not_trimmed.fastq.gz
        echo "No conversion needed" > converted_to_pod5/converted.pod5
    fi
    """
}

process trim_barcodes{
    label 'other_tools'
    publishDir "${params.out_dir}/basecalling_output/", mode: 'copy'
    input:
        path(fastq_not_trimmed) 
    output:
        path("basecalled.fastq.gz"), emit: basecalled_fastq  
    """
    porechop\
     -i basecalled_not_trimmed.fastq.gz\
     -o basecalled.fastq.gz\
     --threads ${params.threads}
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                             Alignment of basecalled sequence to RNA45SN1                                                //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



process align_to_45SN1{
    label 'other_tools'
    publishDir "${params.out_dir}/basecalling_output/", mode:"copy"
    input:
        path(basecalled_fastq) 
        path(reference), stageAs:"reference.fasta"
    output:
        path("filtered.bam"), emit: filtered_bam
        path("filtered.bam.bai"), emit: filtered_bai
        path("filtered.fastq.gz"),emit: filtered_fastq
    """
    minimap2\
     -ax map-ont\
     -t ${params.threads}\
     reference.fasta\
     ${basecalled_fastq}\
     | samtools view -hbS -F 3884\
     | samtools sort\
     > filtered.bam
    #################
    samtools bam2fq filtered.bam --threads ${params.threads} > filtered.fastq
    ##################
    samtools index filtered.bam -@ ${params.threads}
    gzip filtered.fastq -c > filtered.fastq.gz
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                             Filter for read sin fast5 that align 45SN1                                                  //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



process filter_pod5_for_RNA45s_aligning_reads{
    label 'other_tools'
    publishDir "${params.out_dir}/filtered_pod5/", mode:"copy"
    input:
        path(filtered_bam)
        path(filtered_bai)
        path(sample_folder)
        path(converted_pod5), stageAs: "converted.pod5"
    output:
        val 1,emit: done
        path("filtered.pod5"), emit: filtered_pod5
        path("sorted_filtered_reads.txt"), emit: sorted_filtered_read_ids
    """
    mkdir -p filtered_pod5
    (ls ${sample_folder}/*.pod5 && export filetype=pod5) || export filetype=fast5 
    echo \$filetype
    if [ \$filetype == pod5 ]
    then
        python ${projectDir}/bin/filter_pod5.py\
         -i ${filtered_bam}\
         -o .
        pod5 filter ${sample_folder}/*.pod5\
         --ids sorted_filtered_reads.txt\
         --output filtered.pod5\
         --missing-ok\
         --force-overwrite\
         --threads ${params.threads}
    fi
    if [ \$filetype == fast5 ]
    then
        python ${projectDir}/bin/filter_pod5.py\
         -i ${filtered_bam}\
         -o .
        pod5 filter converted.pod5\
         --ids sorted_filtered_reads.txt\
         --output filtered.pod5\
         --missing-ok\
         --force-overwrite\
         --threads ${params.threads}
    fi
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                      Rebasecalling in order to obtain fast5_out option and Event tables in fast5                                        //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process rebasecall_filtered_files{
    label 'dorado_basecaller'
    publishDir "${params.out_dir}/filtered_pod5/", mode:"copy"
    input:
        path(filtered_pod5)
        path(reference), stageAs: "reference.fasta"
        path(basecalling_model)
    output:
        path("filtered_pod5_basecalled.bam") , emit: rebasecalled_bam
        path("filtered_pod5_basecalled.bam.bai"), emit: rebasecalled_bam_bai
    """
    dorado basecaller --estimate-poly-a --emit-moves ${basecalling_model} ${filtered_pod5}\
     | samtools fastq -T "*" --threads ${params.threads}\
     | minimap2 -y --MD -ax map-ont reference.fasta -\
     | samtools sort --threads ${params.threads}\
     | samtools view -b -F 3884 --threads ${params.threads}\
     > filtered_pod5_basecalled.bam 
    samtools index filtered_pod5_basecalled.bam -@ ${params.threads}
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                          Run polyA estimation based on tailfindr                                                        //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process extract_polyA_table{
    label 'other_tools'
    publishDir "${params.out_dir}/taillength_estimation/", mode:"copy"
    input:
        path(rebasecalled_bam)
        path(rebasecalled_bam_bai)
    output:
        val 1,emit: done
        path("tail_estimation.csv")
    """
    python ${projectDir}/bin/extract_polyA_tails.py -i filtered_pod5_basecalled.bam -o .
    """
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                           Run fragment analysis without template to obtain sample specific fragment cluster                             //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



process fragment_analysis_hdbscan{
    label 'other_tools'
    publishDir "${params.out_dir}/fragment_analysis_hdbscan/", mode:"copy"
    input:
        path(filtered_bam)
        path(filtered_bam_bai)
        path(reference), stageAs: "reference.fasta" 
    output:
        path("*")
        val 1, emit: done
    
    """
    python ${projectDir}/bin/fragment_analysis_hdbscan.py -c ${params.threads} -i ${filtered_bam} -r reference.fasta -o ./ -t 0.9 -m 5 -s ${params.color}
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                           Run fragment analysis without template to obtain sample specific fragment cluster                             //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



process fragment_analysis_intensity{
    label 'other_tools'
    publishDir "${params.out_dir}/fragment_analysis_intensity/", mode:"copy"
    input:
        path(filtered_bam)
        path(filtered_bam_bai)
        path(reference), stageAs: "reference.fasta" 
        val(done)
    output:
        path("*")
        val 1, emit: done
    """
    python ${projectDir}/bin/fragment_analysis_intensity.py -c ${params.threads} -i ${filtered_bam} -r reference.fasta -o ./ -t 0.9 -s ${params.color}
    """
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                           Run fragment analysis with templates known to be present in literature                                        //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process template_driven_fragment_analysis{
    label 'other_tools'
    publishDir "${params.out_dir}/template_based_analysis/", mode:"copy"
    input:
        path(filtered_bam)
        path(filtered_bam_bai)
        path(reference), stageAs: "reference.fasta" 
        path(templates), stageAs: "template.csv"
        val(done)
    output:
        path("*")
    """
    python ${projectDir}/bin/template_driven_fragment_analysis.py -c ${params.threads} -i ${filtered_bam} -r reference.fasta -f template.csv -o ./ -t 0.9 -s ${params.color}
    """
}

workflow{
    dorado_basecalling("${params.sample_folder}", "${params.basecalling_model}")
    trim_barcodes(dorado_basecalling.out.fastq_not_trimmed)
    align_to_45SN1(trim_barcodes.out.basecalled_fastq, file("${projectDir}/references/RNA45SN1.fasta"))
    filter_pod5_for_RNA45s_aligning_reads(align_to_45SN1.out.filtered_bam, align_to_45SN1.out.filtered_bai, "${params.sample_folder}", dorado_basecalling.out.converted_pod5)
    rebasecall_filtered_files(filter_pod5_for_RNA45s_aligning_reads.out.filtered_pod5, file("${projectDir}/references/RNA45SN1.fasta"), "${params.basecalling_model}")
    extract_polyA_table(rebasecall_filtered_files.out.rebasecalled_bam, rebasecall_filtered_files.out.rebasecalled_bam_bai)
    fragment_analysis_hdbscan(rebasecall_filtered_files.out.rebasecalled_bam, rebasecall_filtered_files.out.rebasecalled_bam_bai,file("${projectDir}/references/RNA45SN1.fasta"))
    fragment_analysis_intensity(rebasecall_filtered_files.out.rebasecalled_bam, rebasecall_filtered_files.out.rebasecalled_bam_bai,file("${projectDir}/references/RNA45SN1.fasta"),fragment_analysis_hdbscan.out.done)
    template_driven_fragment_analysis(rebasecall_filtered_files.out.rebasecalled_bam, rebasecall_filtered_files.out.rebasecalled_bam_bai,file("${projectDir}/references/RNA45SN1.fasta"), file("${projectDir}/references/Literature_Fragments_and_cut_sites_RNA45SN1.csv"),fragment_analysis_intensity.out.done)
}


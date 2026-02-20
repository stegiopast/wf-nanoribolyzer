#!/usr/bin/ nextflow
//@Grab('com.xlson.groovycsv:groovycsv:1.1')
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
println " "

println "Color"
println params.color
println " "

println "Script folder"
println params.script_folder
println " "

println "Output folder"
println params.out_dir
println " "

println "Basecalling model"
println params.basecalling_model
println " "

println "Threads"
println params.threads
println " "

println "Sample type"
println params.sample_type
println " "

println "Demand"
println params.demand
println " "

println "Model Organism"
println params.model_organism
println " "

println "GPUs in use"
println params.gpus_in_use
println " "

println "Custom reference path"
println params.custom_ref_path
println " "


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
    stageInMode 'symlink'
    input:
    path(sample_folder)
    val basecalling_model

    output:
    val done
    path('basecalling_output/basecalled_not_trimmed.ubam'), emit: basecalled_bam
    path('basecalling_output/sequencing_summary.txt'), optional: true
    path('basecalling_output/basecalled_not_trimmed.fastq.gz'), emit: fastq_not_trimmed
    path('converted_to_pod5/converted.pod5'), emit: converted_pod5, optional: true
    
    script:
    done = 1
    """ 
    mkdir -p basecalling_output
    mkdir -p converted_to_pod5
    
    if compgen -G "${sample_folder}/*.pod5" > /dev/null; then
        export filetype=pod5
    elif compgen -G "${sample_folder}/*.fast5" > /dev/null; then
        export filetype=fast5
    else
        export filetype=bam
    fi 

    #########################
    if [ \$filetype == fast5 ]
    then 
        pod5 convert fast5 ${sample_folder}/*.fast5 --output converted_to_pod5/converted.pod5 --force-overwrite
        if [ ${params.sample_type} == "DNA" ]
        then
            dorado basecaller --no-trim --device ${params.gpus_in_use} --estimate-poly-a --emit-moves ${basecalling_model} converted_to_pod5 > basecalling_output/basecalled_not_trimmed.ubam
            samtools fastq -T "*" --threads ${params.threads} basecalling_output/basecalled_not_trimmed.ubam > basecalling_output/basecalled_not_trimmed.fastq
            gzip -c -1 basecalling_output/basecalled_not_trimmed.fastq > basecalling_output/basecalled_not_trimmed.fastq.gz
            dorado summary basecalling_output/basecalled_not_trimmed.ubam > basecalling_output/sequencing_summary.txt
        else
            dorado basecaller --no-trim --device ${params.gpus_in_use} --estimate-poly-a --emit-moves sup,m5C_2OmeC,inosine_m6A_2OmeA,pseU_2OmeU,2OmeG converted_to_pod5 > basecalling_output/basecalled_not_trimmed.ubam
            samtools fastq -T "*" --threads ${params.threads} basecalling_output/basecalled_not_trimmed.ubam > basecalling_output/basecalled_not_trimmed.fastq
            gzip -c -1 basecalling_output/basecalled_not_trimmed.fastq > basecalling_output/basecalled_not_trimmed.fastq.gz
            dorado summary basecalling_output/basecalled_not_trimmed.ubam > basecalling_output/sequencing_summary.txt
        fi
    fi
    ########################
    if [ \$filetype == pod5 ]
    then
        if [ ${params.sample_type} == "DNA" ]
        then
            dorado basecaller --no-trim --device ${params.gpus_in_use} --estimate-poly-a --emit-moves ${basecalling_model} ${sample_folder} > basecalling_output/basecalled_not_trimmed.ubam
            samtools fastq -T "*" --threads ${params.threads} basecalling_output/basecalled_not_trimmed.ubam > basecalling_output/basecalled_not_trimmed.fastq
            gzip -c -1 basecalling_output/basecalled_not_trimmed.fastq > basecalling_output/basecalled_not_trimmed.fastq.gz
            dorado summary basecalling_output/basecalled_not_trimmed.ubam > basecalling_output/sequencing_summary.txt
        else
            dorado basecaller --no-trim --device ${params.gpus_in_use} --estimate-poly-a --emit-moves sup,m5C_2OmeC,inosine_m6A_2OmeA,pseU_2OmeU,2OmeG ${sample_folder} > basecalling_output/basecalled_not_trimmed.ubam
            samtools fastq -T "*" --threads ${params.threads} basecalling_output/basecalled_not_trimmed.ubam > basecalling_output/basecalled_not_trimmed.fastq
            gzip -c -1 basecalling_output/basecalled_not_trimmed.fastq > basecalling_output/basecalled_not_trimmed.fastq.gz
            dorado summary basecalling_output/basecalled_not_trimmed.ubam > basecalling_output/sequencing_summary.txt
        fi
    fi
    #########################
    if [ \$filetype == bam ]
    then
        samtools view --threads ${params.threads} -bh ${sample_folder}/*bam > basecalling_output/basecalled_not_trimmed.ubam
        samtools fastq -T "*" --threads ${params.threads} basecalling_output/basecalled_not_trimmed.ubam | gzip > basecalling_output/basecalled_not_trimmed.fastq.gz
    fi
    """
}

process trim_barcodes{
    label 'other_tools'
    publishDir "${params.out_dir}/basecalling_output/", mode: 'copy'
    stageInMode 'symlink'
    input:
        path(fastq_not_trimmed) 
    output:
        path("basecalled.fastq.gz"), emit: basecalled_fastq  
    script:
    """
    mkdir -p chunks
    zcat ${fastq_not_trimmed} | split -l 8000000 -d -a 4 - chunks/chunk_
    for chunk in chunks/chunk_*; do
        porechop -i "\$chunk" -o "\${chunk}_trimmed" --threads ${params.threads}
        rm \$chunk
        cat "\${chunk}_trimmed" >> basecalled.fastq
        rm "\${chunk}_trimmed"
    done
    gzip basecalled.fastq
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                        Alignment of basecalled sequence to reference seq                                                //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



process align_to_45SN1{
    label 'other_tools'
    publishDir "${params.out_dir}/filtered_pod5/", mode:"copy"
    stageInMode 'symlink'
    input:
        path(basecalled_fastq) 
        path(reference), stageAs:"reference.fasta"
    output:
        path("filtered_pod5_basecalled.bam"), emit: filtered_bam
        path("filtered_pod5_basecalled.bam.bai"), emit: filtered_bam_bai
        path("filtered_pod5_basecalled.fastq.gz"),emit: filtered_fastq
        val 1, emit: done
    """
    minimap2 -y --MD -ax map-ont -t ${params.threads} reference.fasta ${basecalled_fastq} | samtools view -b -F 3884 | samtools sort > filtered_pod5_basecalled.bam
    #################
    samtools fastq -T "*" --threads ${params.threads} filtered_pod5_basecalled.bam > filtered_pod5_basecalled.fastq
    ##################
    samtools index -@ ${params.threads} filtered_pod5_basecalled.bam 
    gzip -c -1 filtered_pod5_basecalled.fastq > filtered_pod5_basecalled.fastq.gz
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
    stageInMode 'symlink'
    input:
        path(filtered_bam)
        path(filtered_bam_bai)
    output:
        val 1,emit: done
        path("tail_estimation.csv"), emit: tail_esimation_csv
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



// process fragment_analysis_hdbscan{
//     label 'other_tools'
//     publishDir "${params.out_dir}/fragment_analysis_hdbscan/", mode:"copy"
//     input:
//         path(filtered_bam)
//         path(filtered_bam_bai)
//         path(reference), stageAs: "reference.fasta" 
//     output:
//         path("alignment_df.csv"), emit: alignment_df
//         path("fragment_df.csv"), emit: fragment_df
//         path("*")
//         val 1, emit: done
    
//     """
//     python ${projectDir}/bin/fragment_analysis_hdbscan.py -c ${params.threads} -i ${filtered_bam} -r reference.fasta -o ./ -t 0.9 -m 5 -s ${params.color} -d ${params.demand}
//     python ${projectDir}/bin/visualize_clustering_performance_intensity_matrix.py -a ./alignment_df.csv -t ./fragment_df_simple.csv -c ${params.color} -o ./
//     """
// }


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
    stageInMode 'symlink'
    input:
        path(filtered_bam)
        path(filtered_bam_bai)
        path(reference), stageAs: "reference.fasta" 
        val(done)
    output:
        path("alignment_df.csv"), emit: alignment_df
        path("fragment_df.csv"), emit: fragment_df
        path("*")
        val 1, emit: done
    """
    python ${projectDir}/bin/fragment_analysis_intensity.py -c ${params.threads} -i ${filtered_bam} -r reference.fasta -o ./ -t 0.9 -s ${params.color} -d ${params.demand}
    python ${projectDir}/bin/visualize_clustering_performance_intensity_matrix.py -a ./alignment_df.csv -t ./fragment_df_simple.csv -c ${params.color} -o ./
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
    stageInMode 'symlink'
    input:
        path(filtered_bam)
        path(filtered_bam_bai)
        path(reference), stageAs: "reference.fasta" 
        path(templates), stageAs: "template.csv"
        val(done)
    output:
        path("template_fragment_df.csv"), emit: template_csv
        path("template_alignment_df.csv"), emit: template_alignment_csv
        path("start_sites_fragment_based.bed"), emit: start_sites
        path("end_sites_fragment_based.bed"), emit: end_sites
        path("*")
        val 1, emit: done
    """
    python ${projectDir}/bin/template_driven_fragment_analysis.py -c ${params.threads} -i ${filtered_bam} -r reference.fasta -f template.csv -o ./ -t 0.9 -s ${params.color} -d ${params.demand}
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                               Run read-tail analysis based on associated templates of template driven fragment analysis                                 //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process fragment_based_readtail_analysis{
    label 'other_tools'
    publishDir "${params.out_dir}/readtail_analysis/", mode:"copy"
    stageInMode 'symlink'
    input:
        path(filtered_bam)
        path(filtered_bam_bai)
        path(template_csv)
        path(fasta_file)
    output:
        path("*")
        val 1, emit: done
    """
    python ${projectDir}/bin/readtail_analysis.py -i ${filtered_bam} -r ${template_csv} -o ./ -f ${fasta_file}
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                      Visualize poly-A taillengths for template base associations                                                        //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_polyA_associated_templates{
    label 'other_tools'
    publishDir "${params.out_dir}/polyA_template_based/", mode:"copy"
    stageInMode 'symlink'
    input:
        path(template_alignment_csv)
        path(tail_esimation_csv)
        path(templates)
        path(reference)
    output:
        path("*")
        val 1, emit: done
    """
    python ${projectDir}/bin/visualize_taillengths.py -a ${template_alignment_csv} -t ${tail_esimation_csv} -c ${params.color} -o ./ -f ${templates} -r ${reference} -m ${params.model_organism}
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                      Visualize poly-A taillengths for instensity based clusters                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_polyA_associated_intensity_clusters{
    label 'other_tools'
    publishDir "${params.out_dir}/polyA_intensity_based_clusters/", mode:"copy"
    stageInMode 'symlink'
    input:
        path(intensity_alignment_csv)
        path(tail_esimation_csv)
        path(templates)
        path(reference)
    output:
        path("*")
        val 1, emit: done
    """
    python ${projectDir}/bin/visualize_taillengths_clustering.py -a ${intensity_alignment_csv} -t ${tail_esimation_csv} -c ${params.color} -o ./ -f ${templates} -r ${reference} -m ${params.model_organism}
    """
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                      Visualize poly-A taillengths for instensity based clusters                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// process visualize_polyA_associated_hdbscan_clusters{
//     label 'other_tools'
//     publishDir "${params.out_dir}/polyA_hdbscan_based_clusters/", mode:"copy"
//     input:
//         path(hdbscan_alignment_csv)
//         path(tail_esimation_csv)
//         path(templates)
//         path(reference)
//     output:
//         path("*")
//         val 1, emit: done
//     """
//     python ${projectDir}/bin/visualize_taillengths_clustering.py -a ${hdbscan_alignment_csv} -t ${tail_esimation_csv} -c ${params.color} -o ./ -f ${templates} -r ${reference} -m ${params.model_organism}
//     """
// }




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                         Visualize and save intensity matrix                                                             //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_intensity_matrix{
    label 'other_tools'
    publishDir "${params.out_dir}/intensity_matrix/", mode:"copy"
    stageInMode 'symlink'
    input:
        path(alignment_df)
        path(templates)
    output:
        path("*")
        val 1, emit: done
    """
    python ${projectDir}/bin/visualize_intensity_matrix.py -a ${alignment_df} -t ${templates} -c ${params.color} -o ./ -m ${params.model_organism}
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                         Visualize modifications (PsU,m6A)                                                               //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_modifications{
    label 'other_tools'
    publishDir "${params.out_dir}/modification_plots/", mode:"copy"
    stageInMode 'symlink'
    input:
        path(filtered_bam)
        path(filtered_bam_bai)
        path(reference)
        path(modifications_bed)
    output:
        path("*")
        val 1, emit: done
    """
    if [ ${params.sample_type} == "RNA" ]
    then
        python ${projectDir}/bin/visualize_modifications.py -b ${filtered_bam} -r ${reference} -m ${modifications_bed} -o ./
    else
        echo "cDNA used" > modifications.log
    fi
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                         Visualize significant cut sites                                                                 //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_cut_sites{
    label 'other_tools'
    publishDir "${params.out_dir}/cut_site_plots/", mode:"copy"
    stageInMode 'symlink'
    input:
        path(start_sites)
        path(end_sites)
        path(reference)
    output:
        path("*")
        val 1, emit: done
    """
    python ${projectDir}/bin/visualize_cut_sites.py -s ${start_sites} -e ${end_sites}  -r ${reference} -o ./
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                         Visualize reference coverage                                                                    //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_reference_coverage{
    label 'other_tools'
    publishDir "${params.out_dir}/coverage_plots/", mode:"copy"
    stageInMode 'symlink'
    input:
        path(template_alignment_csv)
        path(templates)
        path(reference)
    output:
        path("*.png")
        val 1, emit: done
    """
    python ${projectDir}/bin/visualize_rRNA_coverage.py -a ${template_alignment_csv} -f ${templates}  -r ${reference} -o ./ -c ${params.color} -m ${params.model_organism}
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                               Check if all workflows are done                                                           //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process check_all_done {
    label 'other_tools'
    publishDir "${params.out_dir}/report/", mode: "copy"
    stageInMode 'symlink'
    input:
        val visualize_intensity_matrix_done;
        val visualize_polyA_associated_templates;
        val visualize_polyA_associated_intensity_clusters_done;
        // val visualize_polyA_associated_hdbscan_clusters_done;
        val visualize_modifications_done;
        val visualize_cut_sites_done;
        val visualize_reference_coverage_done;
    output:
        val 1, emit: done
    """
    echo "All visualizations done"
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                               Create HTML report                                                                        //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process create_report {
    label 'other_tools'
    publishDir "${params.out_dir}/", mode: "copy"
    stageInMode 'symlink'
    input:
    // Mandatory inputs
    val all_done_verification
    path intensity_matrix_png, stageAs: "intensity_matrix/intensity_matrix.png"
    path intensity_matrix_html, stageAs: "intensity_matrix/intensity_matrix.html"
    path polyA_tails_intermediates_template_html, stageAs: "polyA_template_based/polyA_tails_intermediates_template.html"
    path polyA_tails_intermediates_min_max_html, stageAs: "polyA_template_based/polyA_tails_intermediates_min_max.html"
    path polyA_tails_intermediates_mean_html, stageAs: "polyA_template_based/polyA_tails_intermediates_mean.html"
    path violinplot_taillength_per_intermediate_png, stageAs: "polyA_template_based/violinplot_taillength_per_intermediate.png"
    path cut_sites_html, stageAs: "cut_site_plots/cut_sites.html"
    path polyA_tails_clustering_html, stageAs: "polyA_intensity_based_clusters/polyA_tails_clustering.html"
    path intensity_matrix_intensity_clustering, stageAs: "fragment_analysis_intensity/intensity_matrix.png"
    path violinplot_taillength_per_cluster_png, stageAs: "polyA_intensity_based_clusters/violinplot_taillength_per_cluster.png"
    // path polyA_tails_hdbscan_clustering_html, stageAs: "polyA_hdbscan_based_clusters/polyA_tails_clustering.html"
    // path intensity_matrix_hdbscan_clustering, stageAs: "fragment_analysis_hdbscan/intensity_matrix.png"
    // path violinplot_taillength_per_hdbscan_cluster_png, stageAs: "polyA_hdbscan_based_clusters/violinplot_taillength_per_cluster.png"
    path coverage_plot_general_absolute, stageAs: "coverage_plots/coverage_fragments_absolute.png"
    path coverage_plot_general_relative, stageAs: "coverage_plots/coverage_fragments_relative.png"
    path coverage_plot_fragments_absolute, stageAs:"coverage_plots/coverage_fragments_absolute_all.png"
    path coverage_plot_fragments_absolute, stageAs:"coverage_plots/coverage_total_sample_absolute.png"
    path coverage_plot_fragments_relative, stageAs:"coverage_plots/coverage_total_sample_relative.png"
    // Optional inputs
    path relative_uridine_modification_abundance_html, stageAs: "modification_plots/uridine_modification_abundance.html"
    path relative_adenine_modification_abundance_html, stageAs: "modification_plots/adenine_modification_abundance.html"
    path relative_guanine_modification_abundance_html, stageAs: "modification_plots/guanine_modification_abundance.html"
    path relative_cytosine_modification_abundance_html, stageAs: "modification_plots/cytosine_modification_abundance.html"

    output:
    path("rRNA_report.html")
    
    script:
    """
    python ${projectDir}/bin/html_report.py -d ./ -o ./
    """
}



workflow{
    if (params.model_organism == "Human"){
        fasta_reference_file = "${projectDir}/references/RNA45SN1.fasta"
        ribosomal_intermediates_file = "${projectDir}/references/Literature_Fragments_and_cut_sites_RNA45SN1.csv"
        modification_reference_file = "${projectDir}/references/rRNA_modifications_conv.bed"
        println fasta_reference_file
        println ribosomal_intermediates_file
        println modification_reference_file
    } else if (params.model_organism == "Yeast"){
        fasta_reference_file = "${projectDir}/references/RDN37-1.fa"
        ribosomal_intermediates_file = "${projectDir}/references/Literature_Fragments_and_cut_sites_RDN37-1.csv"
        modification_reference_file = "${projectDir}/references/rRNA_yeast_modifications_conv.bed"
        println fasta_reference_file
        println ribosomal_intermediates_file
        println modification_reference_file
    }
    dorado_basecalling(
        "${params.sample_folder}", 
        "${params.basecalling_model}"
        )
    trim_barcodes(
        dorado_basecalling.out.fastq_not_trimmed
        )
    align_to_45SN1(
        trim_barcodes.out.basecalled_fastq, 
        file(fasta_reference_file)
        )
    extract_polyA_table(
        align_to_45SN1.out.filtered_bam, 
        align_to_45SN1.out.filtered_bam_bai
        )
    // fragment_analysis_hdbscan(
    //     align_to_45SN1.out.filtered_bam, 
    //     align_to_45SN1.out.filtered_bam_bai,
    //     file(fasta_reference_file)
    //     )
    fragment_analysis_intensity(
        align_to_45SN1.out.filtered_bam, 
        align_to_45SN1.out.filtered_bam_bai,
        file(fasta_reference_file),
        // fragment_analysis_hdbscan.out.done
        extract_polyA_table.out.done
        )
    template_driven_fragment_analysis(
        align_to_45SN1.out.filtered_bam, 
        align_to_45SN1.out.filtered_bam_bai,
        file(fasta_reference_file), 
        file(ribosomal_intermediates_file),
        fragment_analysis_intensity.out.done
        )
    fragment_based_readtail_analysis(
        align_to_45SN1.out.filtered_bam,
        align_to_45SN1.out.filtered_bam_bai,
        template_driven_fragment_analysis.out.template_csv,
        file(fasta_reference_file)
        )
    visualize_intensity_matrix(
        template_driven_fragment_analysis.out.template_alignment_csv,
        file(ribosomal_intermediates_file)
        )
    visualize_polyA_associated_templates(
        template_driven_fragment_analysis.out.template_alignment_csv,
        extract_polyA_table.out.tail_esimation_csv,file(ribosomal_intermediates_file),
        file(fasta_reference_file)
        )
    visualize_polyA_associated_intensity_clusters(
        fragment_analysis_intensity.out.fragment_df,
        extract_polyA_table.out.tail_esimation_csv,
        file(ribosomal_intermediates_file),
        file(fasta_reference_file)
        )
    // visualize_polyA_associated_hdbscan_clusters(
    //     fragment_analysis_hdbscan.out.fragment_df,
    //     extract_polyA_table.out.tail_esimation_csv,
    //     file(ribosomal_intermediates_file),
    //     file(fasta_reference_file)
    //     )
    visualize_modifications(
        align_to_45SN1.out.filtered_bam,
        align_to_45SN1.out.filtered_bam_bai,
        file(fasta_reference_file),
        file(modification_reference_file)
        )
    visualize_cut_sites(
        template_driven_fragment_analysis.out.start_sites,
        template_driven_fragment_analysis.out.end_sites,
        file(fasta_reference_file)
        )
    visualize_reference_coverage(
        template_driven_fragment_analysis.out.template_alignment_csv,
        file(ribosomal_intermediates_file),
        file(fasta_reference_file)
    )
    check_all_done(
        visualize_intensity_matrix.out.done,
        visualize_polyA_associated_templates.out.done,
        visualize_polyA_associated_intensity_clusters.out.done,
        // visualize_polyA_associated_hdbscan_clusters.out.done, 
        visualize_modifications.out.done, 
        visualize_cut_sites.out.done,
        visualize_reference_coverage.out.done
        )
    create_report(
            check_all_done.out.done,
            file("${params.out_dir}/intensity_matrix/intensity_matrix.png"),
            file("${params.out_dir}/intensity_matrix/intensity_matrix.html"),
            file("${params.out_dir}/polyA_template_based/polyA_tails_intermediates_template.html"),
            file("${params.out_dir}/polyA_template_based/polyA_tails_intermediates_min_max.html"),
            file("${params.out_dir}/polyA_template_based/polyA_tails_intermediates_mean.html"),
            file("${params.out_dir}/polyA_template_based/violinplot_taillength_per_intermediate.png"),
            file("${params.out_dir}/cut_site_plots/cut_sites.html"),
            file("${params.out_dir}/polyA_intensity_based_clusters/polyA_tails_clustering.html"),
            file("${params.out_dir}/fragment_analysis_intensity/intensity_matrix.png"),
            file("${params.out_dir}/polyA_intensity_based_clusters/violinplot_taillength_per_cluster.png"),
            // file("${params.out_dir}/polyA_hdbscan_based_clusters/polyA_tails_clustering.html"),
            // file("${params.out_dir}/fragment_analysis_hdbscan/intensity_matrix.png"),
            // file("${params.out_dir}/polyA_hdbscan_based_clusters/violinplot_taillength_per_cluster.png"),
            file("${params.out_dir}/coverage_plots/coverage_fragments_absolute.png"),
            file("${params.out_dir}/coverage_plots/coverage_fragments_relative.png"),
            file("${params.out_dir}/coverage_plots/coverage_fragments_absolute_all.png"),
            file("${params.out_dir}/coverage_plots/coverage_total_sample_absolute.png"),
            file("${params.out_dir}/coverage_plots/coverage_total_sample_relative.png"),
            file("${params.out_dir}/modification_plots/uridine_modification_abundance.html"),
            file("${params.out_dir}/modification_plots/adenine_modification_abundance.html"),
            file("${params.out_dir}/modification_plots/guanine_modification_abundance.html"),
            file("${params.out_dir}/modification_plots/cytosine_modification_abundance.html")
    )
}



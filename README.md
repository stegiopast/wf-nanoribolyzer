# NanoRibolyzer

[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md) ![GitHub Release](https://img.shields.io/github/v/release/dietvin/pod5Viewer) ![PyPI - Version](https://img.shields.io/pypi/v/pod5Viewer) ![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/dietvin/pod5Viewer/total)

## Abstract
NanoRibolyzer is an Epi2Me compatible tool to analyse ribosomal RNA biogenesis pathway at single nucleotide resolution. 
Nanoribolyzer aligns reads to the human RNA45SN1 of hg38. A template based read association approach allows for quantification of known ribosomal intermediates.
Several template free clustering approaches can are integrated to detect and study unknown ribosomal RNA intermediates. 
In addition to the template association NanoRibolyzer performs polyA tail estimation, finds abundant cut sites, extracts 5' terminal base sequences for motif analysis 
and detects the relativa abundance of specific RNA modifications. Outputs of NanoRibolyzer can be used to disect specific subpopulations of ribosomal RNA reads and asses characteristic properties   
on a single nucleotide resolution level. 

## Userguide
To use NanoRibolyzer [Epi2Me](https://labs.epi2me.io/downloads/),[Nextflow](https://www.nextflow.io/docs/latest/install.html) and [Docker](https://docs.docker.com/engine/install/) must be installed.

To download NanoRibolyzer open Epi2Me and navigate to Launch and press the button Import Workflow. A pop-up window will appear in which you should copy the following link: "https://github.com/stegiopast/wf-nanoribolyzer"
Press download and the Workflow should be integrated in Epi2Me. On Windows we recommend to install [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) and subsequently install Epi2Me, Nextflow and Docker within the command line.
Also make sure to install the [nvidia-container-toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) on WSL.

## Methods

![General Pipeline](./figures/General_pipeline.png)

The analysis workflow of NanoRibolyzer starts with the pod5 output format of ONT’s MinKnow. Reads are basecalled using dorado basecaller (cite). All sequenced reads are basecalled and trimmed with Porechop (cite). The trimmed reads become aligned with the map-ont flag of minimap2 (cite) to the 45SN1 reference of hg38 (cite). The ids of 45SN1 aligning reads are used to filter the original pod5 file. The filtered pod5 file is rebasecalled using the integrated models for modification detection and polyA taillength of dorado. The read ids in the resulting unaligned bam file is used to collect metainformation about reads on a single nucleotide resolution. 
Rebasecalled reads become aligned to the 45SN1 reference. Resulting bam files are used to perform several clustering algorithms. The pipeline includes a reference-based and two reference-free clustering approaches.  

# Template based fragment association

![Template based clustering](./figures/Template_based_clustering.png)

The reference-based algorithm performs an association of sequenced reads to literature based ribosomal intermediates. In this approach the minimal overlap pairs between a query read and all possible intermediates becomes determined. The minimal overlap is determined by defining the minimal relative overlap of query over intermediate and intermediate over query. After defining the minimal overlap for each possible query-intermediate pair, the argument of the maximal overlap-pair is used for the intermediate association. Read clusters are stored in tsv format including the read ids, the absolute and relative amount and the start and end sites of all reads in a cluster. Non clustered reads are stored in a separate file. Moreover, bed files with the read clusters are stored to enable visualization on the integrate genome viewer (igv). The 45SN1 reference fasta in the ribolyzer references repository must be used for visualization. 


# Template free fragment association

![Template free clustering](./figures/Template_free_clustering.png)

The reference-free algorithms are based on a preceding intensity matrix construction. Read ids become embedded in a 45SN1 length x 45SN1 length matrix by using their start and end points of the alignment as coordinates. The intensity of a coordinate in the matrix is determined by the number of reads aligning to it. 

1. In the first clustering approach all reads sharing a start and end site are interpreted as a read cluster. The approach is most performant with reads not underlying degradation processes. 
2. The second approach is using a hierarchical density-based clustering approach with alternating neighbourhood (HDBSCAN) to cluster reads sharing similar start and end sites in the intensity matrix. HDBSCAN determines clusters of high intensity read groups having several neighbouring read groups in the intensity matrix. Coordinates on the intensity matrix with a high intensity lacking neighbours are defined as independent clusters. Resulting clusters of the HDBSCAN approach can be summarized by either constructing a consensus sequence of reads belonging to a cluster (higher demand) or by extracting a reference sequence from the 45SN1 of hg38 by using the minimal start point and the maximal end point of the cluster. (lower demand)  

# Poly-A estimation
NanoRibolyzer is determining the polyA lengths with the integrated polyA length estimation tool of dorado. Due to the experimental design of the Oxford Nanopore Technologies related sequencing library, polyA tails can only be captured in the entire length when performing directRNA sequencing approaches, since these protocols exclusively include a ligation of oligo_dT primers at the 3'end of RNA fragments. NanoRibolyzer performs polyA estimation on cDNA as well, but it is crucial to consider the oligo_dT primers aligning on any sterical possible position of polyA tail. 

# Modification detection
Modification detection is performed with the integrated modification detection models of dorado. Modifications can only be detected when performing directRNA sequencing approaches, since reverse transcription and strand switching to obtain cDNA erases the chemical signature of RNA modifications.    

# Cut-site determination


![Cut site determination](./figures/Determine_cut_sites.png)

A determination of significant abundant start and end sites is computed using the start and end site information of bam files. Start and end sites with absolute abundance are determined over the 45SN1 reference. For each reference-based intermediate cut site interval, the mean and standard-deviation (stdd) of start and end sites become determined respectively. Cut sites occurring at least mean + (2*stdd) times are considered as significant abundant cut-sites. For overlapping reference-based intermediate cut site intervals, metrics are determined by the mean of means and the mean of stdds. Output files are stored in tsv and bed-file format including information about their relative abundance and the cut site location. The bed files can be visualized using igv.  

# Extraction of 20 nucleotides before polyA tail

![Motif extraction](./figures/Motif_extraction.png)

NanoRibolyzer performs an extraction of the last 20 nucleotides before a polyA tail. Each read is respectively investigated. Starting from the 3’ end of a trimmed read the polyA tail starting position is determined by finding the first triple adenine (AAA) pattern on the read. Starting from the determined point the algorithm iterates over the nucleotides towards the 5’end until less or equal 50% of the nucleotides are adenines (A). From that point the algorithm iterates towards the 3’end again, until the observed position of the iteration is an adenine and more or equal 90% of the nucleotides are adenines to define the end point of the polyA tail. The last 20 nucleotides from 5’ to 3’ end of the determined end point are extracted and stored in three ways. First the end sites are sorted in a txt format file including only the 20-nucleotide sequence, which can be used as direct input for motif [enrichment analysis tools](https://meme-suite.org/meme/tools/streme). Secondly, the reads are stored with their read ids and the 20-nucleotide sequence in fasta format. Thirdly, the reads are stored with their read ids, the 20-nucleotide sequence and the 20-nucleotide sequence + trashed 3’end sequence for extraction assessment in a second fasta file.

## Outputs

NanoRibolyzer provides files, tables and graphics in different subfolders. Here we describe all the outputs in detail to enable users the performance of downstream analysis. 
All the outputs will be provided in the default workfolder of Epi2Me.

```bash
├── basecalling_output      # Data for basecalling, trimming and first alignment to 45SN1 of hg38
│   │
│   ├── basecalled_not_trimmed.bam              # Unaligned basecalling output of dorado untrimmed
│   ├── basecalled_not_trimmed.fastq.gz         # Converted fastq file of basecalling output untrimmed 
│   ├── basecalled.fastq.gz                     # Converted fastq file of basecalling output trimmed with porechop
│   ├── filtered.bam                            # Alignment of basecalled.fastq.gz to 45SN1 of hg38 
│   ├── filtered.bam.bai                        # Index for above
│   ├── filtered.fastq.gz                       # Converted fastq file of filtered.bam 
│   └── sequencing_summary.txt                  # Summary fuile of dorado basecalling
│
├── converted_to_pod5                           # If original seqeuncing files were written in fast5 format they become converted into pod5 first and will be stored in this directory
│   │
│   └── converted.pod5                          # File is written only if original sequencing file was in fast5 format
│
├── filtered_pod5
│   │
│   ├── filtered.pod5                           # Pod5 file with files aligning to 45SN1 of hg38
│   ├── filtered_pod5_basecalled.bam            # 45SN1 aligned basecalling output of dorado from filtered.pod5 untrimmed. This file includes move table (all) and modification tags (directRNA). 
│   ├── filtered_pod5_basecalled.bam.bai        # Index for above
│   └── sorted_filtered_reads.txt               # List with read_ids aligning 45SN1 of hg38
│
├── intensity_matrix                            # Folder contains
│   │
│   ├── intensity_matrix.csv                    # Table with read ids and number of reads aligning to (start site, end site) pairs  
│   ├── intensity_matrix.html                   # Interactive html file with 50000 most abundant (start site, end site) pairs
│   └── intensity_matrix.png                    # An image of the intensity matrix used for the clustering approach
│
├── template_based_analysis
│   │
│   ├── template_alignment_df.csv               # List of single aligned reads with aligned sequences reconstructed by cigarstrings, includes start and end sites on the reference. The percentage of overlap with associated fragment is listed. 
│   ├── template_fragment_df.csv                # List with read_ids being associated to fragments from literature. Table stores the extracted sequence of the literature fragment and the start and end site on 45SN1.  
│   ├── cutting_sites_general.csv               # Cut sites occuring in significatn abundance in relation to the mean abundance fo reads on 45SN1
│   ├── cutting_sites_fragment_based.csv        # Cut sites occuring significant abundance in relation to mean abundance of the literature fragments they are associated to. Threshholds for overlapping fragments were determined by the mean of means. 
│   ├── start_sites_general.bed                 # Bed file for visualization of general start sites in IGV
│   ├── end_sites_general.bed                   # Bed file for visualization of general end sites in IGV
│   ├── start_sites_fragment_based.bed          # Bed file for visualization of fragment based start sites in IGV
│   ├── end_sites_fragment_based.bed            # Bed file for visualization of fragment based end sites in IGV
│   ├── template_driven.bed                     # Bed file visualizing fragment abundance in igv
│   ├── most_abundant_template_driven.bed       # Bed file visualizing fragment abundance in igv. Only fragments carrying more than mean * 2.stdd reads.
│   └── template_driven_analysis.bed            # Log file of template_based_analysis script
│
├── fragment_hdbscan_analysis                   # Outputs of template free clustering approaches. Same structure is given for fragment_intensity_analysis.
│   │
│   ├── alignment_df.csv                        # List of aligned reads with aligned sequences reconstructed by cigarstrings, includes start and end sites on the reference
│   ├── fragment_df.csv                         # For each fragment being determined by template free association a consensus sequence, proportional sequence, read ids, absolute and relative number of reads are provided
│   ├── fragment_df_simple.csv                  # A simpler form of the table above
│   ├── unclustered_reads_df.csv                # Stores all reads, which could not be clustered in a read group. The table has an equal format as alignment_df.csv  
│   ├── no_template.bed                         # Representation of all fragments found with template free clustering approaches.    
│   ├── most_abundant_no_template.bed           # Fragments found with template free clustering approaches carrying more than mean * 2stdd of the reads  
│   ├── top_150_no_template.bed                 # Top 150 Fragments found with template free clustering carrying most of the reads
│   ├── intensity_matrix.png                    # Intensity matrix used for the clustering approach with 300 clusters carrying the most reads. Clusters are shown as red rectangles of which the upper left corner represents common (start,end) sites. 
│   └── fragment_analysis.log                   # Log file of fragment_hdbscan_analysis script
│
├── readtail_analysis                           # Extract the last 20 nucleotides in every reads before the polyA tail. Output of this folder can be used for motif enrichment analysis
│   │
│   ├── X.txt                                   # List of 20 nucleotide long sequences of readtails before polyA tail for all reads associated with literature fragment (X).  
│   ├── X.fasta                                 # Fasta of 20 nucleotide long sequences of readtails before polyA tail for all reads associated with literature fragment (X). Header indicates the read id.
│   ├── X_tail_comparison.fasta                 # Fasta of 20 nucleotide long sequences of readtails before polyA tail for all reads associated with literature fragment (X). Header indicates the read id. Additionally shows the excluded polyA tail.
│   └── X.json                                  # Dictionairy with modet abundant nucleotides for each position of the last 20 nucleotides before polyA tail for all reads associated with literature fragment (X).
│
├── taillength_estimation                       # Output folder of polyA taillength assessment
│   │
│   ├── tail_estimation.csv                     # Lists taillength for each read_id aligned to 45SN1
│
├── polyA_template_based                        # Assesment of polyA taillengths for reads associated to specific literature based fragments
│   │                  
│   ├── polyA_tails_intermediates_mean.html     # PolyA taillength of reads associated with literature based fragments. Start and end sites of shown fragments are determined by mean start and end sites of reads associated to a fragment. 
│   ├── polyA_tails_intermediates_min_max.html  # PolyA taillength of reads associated with literature based fragments. Start and end sites of shown fragments are determined by minimal start and maximal end site of reads associated to a fragment.      
│   ├── polyA_tails_intermediates_template.html # PolyA taillength of reads associated with literature based fragments. Start and end sites fit to the start and end sites of literature based fragments.
│   ├── taillength_per_intermediate_mean.csv    # Table storing the information for the html objects above
│   ├── taillength_per_intermediate_min_max.csv # Table storing the information for the html objects above
│   └── violinplot_taillength_per_intermediate.png   # Figure showing the polyA taillength for all reads associated to specific literature based fragments 
│
├── polyA_hdbscan_based_clusters                # Assesment of polyA taillengths for reads associated to specific fragments resulting from hdbscan. Same folder structure for polyA_intensity_base_clusters
│   │                   
│   ├── polyA_tails_clustering.html             # PolyA taillength of reads associated with reference free fragments. The 300 most abundant fragments are shown.  
│   ├── polyA_tails_intermediates_min_max.html  # PolyA taillength of reads associated with literature based fragments. Start and end sites of shown fragments are determined by minimal start and maximal end site of reads associated to a fragment.      
│   └── polyA_tails_intermediates_template.html # PolyA taillength of reads associated with literature based fragments. Start and end sites fit to the start and end sites of literature based fragments.
│
├── coverage_plots                              # Plots showing the read coverage of reference based templated after association
│   │
│   ├── coverage_fragments_absolute.png         # Coverage plot for each reference based fragment after read associtaion. Normalize to the maximum expression of a respective fragment.
│   ├── coverage_fragments_absolute_all.png     # Coverage plot for each reference based fragment after read associtaion. Normalize to the number of all reads in the sample. 
│   ├── coverage_fragments_relative.png         # Coverage plot for each reference based fragment after read associtaion. Normalization to the maximum expression of a respective fragment and shown as percentage.
│   ├── coverage_total_sample_absolute.png      # Read coverage of the whole 45SN1 normalized to number of all reads in the sample.
│   └── coverage_total_sample_relative.png      # Read coverage of the whole 45SN1 normalized to number of all reads in the sample and shown as percentage.
│
├── cut_site_plot                               # Visualization for significantly abundant cut sites
│   │
│   └── cut_sites.html                          # Cut sites determined by fragment based threshold cut site determination of template based analysis. Normalized to the number of all reads in the sample.
│
└── rRNA_report.html                            # HTML report including all plots decribed above. Will be automatically integrated in Epi2Me.
```

## Software versions
```bash
name: other_tools
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - _libgcc_mutex=0.1
  - _openmp_mutex=4.5
  - asttokens=2.4.1
  - attrs=23.1.0
  - aws-c-auth=0.7.0
  - aws-c-cal=0.6.0
  - aws-c-common=0.8.23
  - aws-c-compression=0.2.17
  - aws-c-event-stream=0.3.1
  - aws-c-http=0.7.11
  - aws-c-io=0.13.28
  - aws-c-mqtt=0.9.0
  - aws-c-s3=0.3.13
  - aws-c-sdkutils=0.1.11
  - aws-checksums=0.1.16
  - aws-crt-cpp=0.21.0
  - aws-sdk-cpp=1.10.57
  - backcall=0.2.0
  - backports=1.1
  - backports.functools_lru_cache=1.6.5
  - bcftools=1.20
  - beautifulsoup4=4.12.3
  - bedtools=2.31.1
  - biopython=1.78
  - blas=1.0
  - bokeh=3.4.1
  - bottleneck=1.3.7
  - brotli=1.0.9
  - brotli-bin=1.0.9
  - brotli-python=1.0.9
  - bzip2=1.0.8
  - c-ares=1.19.1
  - ca-certificates=2024.6.2
  - cached-property=1.5.2
  - cached_property=1.5.2
  - certifi=2024.6.2
  - charset-normalizer=3.2.0
  - chart-studio=1.1.0
  - click=8.1.7
  - cloudpickle=3.0.0
  - comm=0.2.1
  - contourpy=1.2.0
  - coreutils=8.32
  - curl=8.7.1
  - cycler=0.11.0
  - cyrus-sasl=2.1.28
  - cytoolz=0.12.2
  - dask=2024.5.0
  - dask-core=2024.5.0
  - dask-expr=1.1.0
  - dbus=1.13.18
  - debugpy=1.8.0
  - decorator=5.1.1
  - distributed=2024.5.0
  - exceptiongroup=1.2.0
  - executing=2.0.1
  - expat=2.6.2
  - font-ttf-dejavu-sans-mono=2.37
  - font-ttf-inconsolata=2.001
  - font-ttf-source-code-pro=2.030
  - font-ttf-ubuntu=0.83
  - fontconfig=2.14.1
  - fonts-anaconda=1
  - fonts-conda-forge=1
  - fonttools=4.51.0
  - freetype=2.12.1
  - fsspec=2024.3.1
  - gdbm=1.18
  - gflags=2.2.2
  - giflib=5.2.1
  - glib=2.78.4
  - glib-tools=2.78.4
  - glog=0.6.0
  - gsl=2.7.1
  - gst-plugins-base=1.14.1
  - gstreamer=1.14.1
  - gxf2bed=0.2.2
  - hdf5=1.12.1
  - heapdict=1.0.1
  - htslib=1.20
  - icu=73.1
  - idna=3.7
  - importlib-metadata=7.0.1
  - importlib_metadata=7.0.1
  - importlib_resources=6.1.1
  - intel-openmp=2023.1.0
  - ipykernel=6.28.0
  - ipython=8.16.1
  - jedi=0.19.1
  - jinja2=3.1.4
  - joblib=1.4.2
  - jpeg=9e
  - jsonschema=4.19.2
  - jsonschema-specifications=2023.7.1
  - jupyter_client=8.6.0
  - jupyter_core=5.5.0
  - k8=0.2.5
  - kaleido-core=0.2.1
  - keyutils=1.6.1
  - kiwisolver=1.4.4
  - kneed=0.8.5
  - krb5=1.20.1
  - lcms2=2.12
  - ld_impl_linux-64=2.40
  - lerc=3.0
  - libabseil=20230125.3
  - libaec=1.0.6
  - libarrow=12.0.1
  - libblas=3.9.0
  - libbrotlicommon=1.0.9
  - libbrotlidec=1.0.9
  - libbrotlienc=1.0.9
  - libcblas=3.9.0
  - libclang=14.0.6
  - libclang13=14.0.6
  - libcrc32c=1.1.2
  - libcups=2.4.2
  - libcurl=8.7.1
  - libdeflate=1.20
  - libedit=3.1.20230828
  - libev=4.33
  - libevent=2.1.12
  - libffi=3.4.4
  - libgcc=7.2.0
  - libgcc-ng=13.2.0
  - libgfortran=3.0.0
  - libgfortran-ng=13.1.0
  - libgfortran5=13.1.0
  - libglib=2.78.4
  - libgomp=13.2.0
  - libgoogle-cloud=2.12.0
  - libgrpc=1.56.2
  - libhwloc=2.9.1
  - libiconv=1.17
  - liblapack=3.9.0
  - libllvm14=14.0.6
  - libnghttp2=1.57.0
  - libnl=3.9.0
  - libnuma=2.0.18
  - libopenblas=0.3.21
  - libpng=1.6.39
  - libpq=12.17
  - libprotobuf=4.23.3
  - libsodium=1.0.18
  - libsqlite=3.46.0
  - libssh2=1.11.0
  - libstdcxx-ng=13.2.0
  - libthrift=0.18.1
  - libtiff=4.2.0
  - libutf8proc=2.8.0
  - libuuid=1.41.5
  - libwebp=1.3.2
  - libwebp-base=1.3.2
  - libxcb=1.15
  - libxkbcommon=1.0.1
  - libxml2=2.10.4
  - libzlib=1.2.13
  - locket=1.0.0
  - lz4=4.3.2
  - lz4-c=1.9.4
  - mappy=2.28
  - markupsafe=2.1.3
  - mathjax=2.7.5
  - matplotlib=3.8.4
  - matplotlib-base=3.8.4
  - matplotlib-inline=0.1.6
  - minimap2=2.28
  - mkl=2023.1.0
  - mkl-service=2.4.0
  - mkl_fft=1.3.8
  - mkl_random=1.2.4
  - msgpack-python=1.0.7
  - mysql=5.7.24
  - nbformat=5.9.2
  - ncurses=6.4
  - nest-asyncio=1.6.0
  - nextflow=24.04.2
  - nspr=4.35
  - nss=3.89.1
  - numexpr=2.8.7
  - numpy-base=1.26.4
  - ont-fast5-api=4.1.3
  - openblas=0.3.4
  - openjdk=11.0.13
  - openjpeg=2.4.0
  - openssl=3.3.1
  - orc=1.9.0
  - packaging=23.2
  - pandas=2.2.2
  - parso=0.8.3
  - partd=1.4.1
  - pcre=8.45
  - pcre2=10.42
  - perl=5.34.0
  - pexpect=4.8.0
  - pickleshare=0.7.5
  - pillow=10.3.0
  - pip=24.0
  - platformdirs=3.10.0
  - plotly=5.22.0
  - ply=3.11
  - pooch=1.7.0
  - porechop=0.2.4
  - prompt-toolkit=3.0.43
  - prompt_toolkit=3.0.43
  - psutil=5.9.5
  - ptyprocess=0.7.0
  - pure_eval=0.2.2
  - pyarrow-hotfix=0.6
  - pybedtools=0.10.0
  - pybind11-abi=5
  - pygments=2.16.1
  - pyparsing=3.0.9
  - pyqt=5.15.10
  - pyqt5-sip=12.13.0
  - pysam=0.22.1
  - pysocks=1.7.1
  - python=3.9.19
  - python-dateutil=2.9.0post0
  - python-fastjsonschema=2.16.2
  - python-lmdb=1.4.1
  - python-tzdata=2023.3
  - python_abi=3.9
  - pytz=2024.1
  - pyyaml=6.0.1
  - pyzmq=25.1.2
  - qt-main=5.15.2
  - rdma-core=51.1
  - re2=2023.03.02
  - readline=8.2
  - referencing=0.30.2
  - requests=2.32.2
  - retrying=1.3.3
  - rpds-py=0.10.6
  - s2n=1.3.46
  - samtools=1.19.2
  - scikit-learn=1.4.2
  - scipy=1.12.0
  - seaborn=0.13.2
  - seqkit=2.8.2
  - sip=6.7.12
  - six=1.16.0
  - snappy=1.1.10
  - sortedcontainers=2.4.0
  - sqlite=3.45.3
  - stack_data=0.6.2
  - tar=1.34
  - tbb=2021.8.0
  - tblib=2.0.0
  - threadpoolctl=3.1.0
  - tk=8.6.14
  - tomli=2.0.1
  - toolz=0.12.0
  - tornado=6.3.3
  - traitlets=5.14.3
  - typing-extensions=4.11.0
  - typing_extensions=4.11.0
  - tzdata=2024a
  - ucx=1.14.1
  - unicodedata2=15.1.0
  - urllib3=2.2.1
  - vsearch=2.28.1
  - wcwidth=0.2.8
  - wheel=0.43.0
  - xyzservices=2023.10.1
  - xz=5.4.6
  - yaml=0.2.5
  - zeromq=4.3.5
  - zict=3.0.0
  - zipp=3.17.0
  - zlib=1.2.13
  - zstd=1.5.5
  - pip:
      - fastq==2.0.4
      - h5py==3.10.0
      - iso8601==2.0.0
      - kaleido==0.2.1
      - lib-pod5==0.3.10
      - more-itertools==9.1.0
      - numpy==1.25.1
      - pod5==0.3.10
      - polars==0.20.31
      - progressbar33==2.4
      - pyarrow==14.0.2
      - pyside6==6.7.1
      - pyside6-addons==6.7.1
      - pyside6-essentials==6.7.1
      - setuptools==70.0.0
      - shiboken6==6.7.1
      - tenacity==8.2.3
      - tqdm==4.65.0
      - vbz-h5py-plugin==1.0.1
prefix: /home/stefan/anaconda3/envs/other_tools
```

## References

## Funding


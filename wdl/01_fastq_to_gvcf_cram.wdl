workflow sentieon_ccdg_fastq_vcf {
  # Inputs
  Array [String] file_fastq_r1s
  Array [String] file_fastq_r2s

  Array[String] read_groups
  String sample_name

  # Known sites
  File dbsnp_vcf = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
  File dbsnp_index = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
  File ref_alt = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt"
  File ref_fasta = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
  File ref_fasta_index = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
  File ref_dict = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict"
  File ref_bwt = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt"
  File ref_sa = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa"
  File ref_amb = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb"
  File ref_ann = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann"
  File ref_pac = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac"
  File mills_vcf =  "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  File mills_vcf_index  =  "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
  File hapmap_vcf =  "gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz"
  File hapmap_vcf_index =  "gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi"
  File omni_vcf =  "gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz"
  File omni_vcf_index =  "gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi"
  File onekg_vcf =  "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  File onekg_vcf_index =  "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
  File onekg_indel_vcf =  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
  File onekg_indel_vcf_index =  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
  File axiom_poly_vcf = "gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
  File axiom_poly_vcf_index = "gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"
 
  Array[File] bqsr_vcfs = [dbsnp_vcf,mills_vcf,onekg_indel_vcf] 
  Array[File] bqsr_tbis = [dbsnp_index,mills_vcf_index,onekg_indel_vcf_index]


  # Workflow configurations
  ## Optional workflow stages
  Boolean output_align_file = true # If false, output cram will be empty
  Boolean upload_gvcf = false
  Boolean output_gvcf = true
  Boolean output_vcf = true
  ## BQSR intervals
  String bqsr_intervals = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
  ## Readwriter read_filter args, starting after the QualCalFilter table
  String readwriter_readfilter_args = ",prior=-1.0,indel=false,levels=10/20/30,min_qual=6"
  ## Variant calling algorithm
  String calling_algo = "Haplotyper"
  String calling_intervals = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
  ## Extra driver parameters
  String lc_driver_args = "--traverse_param=200000/10000"
  String dedup_driver_args = "--traverse_param=200000/10000"
  String readwriter_driver_args = ""
  String bqsr_driver_args = ""
  String calling_driver_args = ""
  String genotyping_driver_args = ""
  ## Extra algo parameters
  String bwa_args = "-Y -h 0,0"
  String bwa_chunk_size = "10000000"
  String sort_args = "--block_size 512M --bam_compression 1"
  String lc_args = ""
  String dedup_args = ""
  String readwriter_args = "--cram_write_options version=3.0,compressor=rans"
  String bqsr_args = ""
  String calling_args = ""
  String genotyping_args = ""
  ## Alignment file formats
  Boolean alignment_cram = false
  Boolean dedup_cram = false
  Boolean readwriter_cram = true

  # Sentieon License configuration
  String sentieon_license_server = ""
  Boolean use_instance_metadata = false
  String sentieon_auth_mech = ""
  String sentieon_license_key = ""

  # Execution configuration
  String disk
  String threads = "64"
  String memory = "55 GB"
  Int preemptible_tries = 2
  String sentieon_version = "201808.06"
  String docker = "dnastack/sentieon-google-cloud:${sentieon_version}"
  String sentieon_release_dir = "/opt/sentieon/sentieon-genomics-${sentieon_version}"
  
  call SentieonFastqToVcf {
    input:
      # Inputs
      file_fastq_r1s = file_fastq_r1s,
      file_fastq_r2s = file_fastq_r2s,
      read_groups = read_groups,
      sample_name = sample_name,
      # Known sites
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_index = dbsnp_index,
      bqsr_vcfs = bqsr_vcfs,
      bqsr_tbis = bqsr_tbis,
      mills_vcf = mills_vcf,
      mills_vcf_index = mills_vcf_index,
      hapmap_vcf = hapmap_vcf,
      hapmap_vcf_index = hapmap_vcf_index,
      omni_vcf = omni_vcf,
      omni_vcf_index = omni_vcf_index,
      onekg_vcf = onekg_vcf,
      onekg_vcf_index = onekg_vcf_index,
      axiom_poly_vcf = axiom_poly_vcf,
      axiom_poly_vcf_index = axiom_poly_vcf_index,
      # Reference files
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_alt = ref_alt,
      ref_bwt = ref_bwt,
      ref_sa = ref_sa,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_pac = ref_pac,
      # Workflow configurations
      ## Optional workflow stages
      output_align_file = output_align_file,
      upload_gvcf = upload_gvcf,
      output_gvcf = output_gvcf,
      output_vcf = output_vcf,
      ## BQSR intervals
      bqsr_intervals = bqsr_intervals,
      ## BQSR read_filter args, starting with QualCalFilter
      readwriter_readfilter_args = readwriter_readfilter_args,
      ## Variant calling parameters
      calling_algo = calling_algo,
      calling_intervals = calling_intervals,
      ## Extra driver parameters
      lc_driver_args = lc_driver_args,
      dedup_driver_args = dedup_driver_args,
      readwriter_driver_args = readwriter_driver_args,
      bqsr_driver_args = bqsr_driver_args,
      calling_driver_args = calling_driver_args,
      genotyping_driver_args = genotyping_driver_args,
      ## Extra algo parameters
      bwa_args = bwa_args,
      bwa_chunk_size = bwa_chunk_size,
      sort_args = sort_args,
      lc_args = lc_args,
      dedup_args = dedup_args,
      readwriter_args = readwriter_args,
      bqsr_args = bqsr_args,
      calling_args = calling_args,
      genotyping_args = genotyping_args,
      ## Alignment file formats
      alignment_cram = alignment_cram,
      dedup_cram = dedup_cram,
      readwriter_cram = readwriter_cram,
      # Sentieon License configuration
      sentieon_license_server = sentieon_license_server,
      use_instance_metadata = use_instance_metadata,
      sentieon_auth_mech = sentieon_auth_mech,
      sentieon_license_key = sentieon_license_key,
      # Execution configuration
      disk = disk,
      threads = threads,
      memory = memory,
      preemptible_tries = preemptible_tries,
      docker = docker,
      sentieon_release_dir = sentieon_release_dir
  }

  call stripTags {
    input:
      cram = SentieonFastqToVcf.alignment,
      cram_index = SentieonFastqToVcf.alignment_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      sample_name = sample_name,
      docker = docker
  }

  call cramStats {
    input:
      cram = SentieonFastqToVcf.alignment,
      cram_index = SentieonFastqToVcf.alignment_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      sample_name = sample_name,
      docker = docker
  }

  output {
    File preprocessed_alignment = stripTags.stripped_cram
    File preprocessed_alignment_index = stripTags.stripped_cram_index
    File gvcf = SentieonFastqToVcf.gvcf
    File gvcf_index = SentieonFastqToVcf.gvcf_index
    File dedup_metrics = SentieonFastqToVcf.dedup_metrics
    File mq_metrics = SentieonFastqToVcf.mq_metrics
    File mq_plot = SentieonFastqToVcf.mq_plot
    File qd_metrics = SentieonFastqToVcf.qd_metrics
    File qd_plot = SentieonFastqToVcf.qd_plot
    File gc_summary = SentieonFastqToVcf.gc_summary
    File gc_metrics = SentieonFastqToVcf.gc_metrics
    File gc_plot = SentieonFastqToVcf.gc_plot
    File as_metrics = SentieonFastqToVcf.as_metrics
    File is_metrics = SentieonFastqToVcf.is_metrics
    File is_plot = SentieonFastqToVcf.is_plot
    File stats = cramStats.stats
  }
  
  meta {
    author: "Heather Ward"
    email: "heather@dnastack.com"
    description: "Data analysis pipeline used for converting fastqs to cram and gvcf for the MSSNG DB6 release"
  }
}


task SentieonFastqToVcf {
  # Inputs
  Array[File] file_fastq_r1s
  Array[File] file_fastq_r2s
  Array[String] read_groups
  String sample_name

  # Known sites
  File? dbsnp_vcf
  File? dbsnp_index
  Array[File] bqsr_vcfs
  Array[File] bqsr_tbis
  
  File mills_vcf
  File mills_vcf_index
  File hapmap_vcf
  File hapmap_vcf_index
  File omni_vcf
  File omni_vcf_index
  File onekg_vcf
  File onekg_vcf_index
  File axiom_poly_vcf
  File axiom_poly_vcf_index

  # Reference files
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_alt
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

  # Workflow configurations
  ## Optional workflow stages
  Boolean output_align_file
  Boolean upload_gvcf
  Boolean output_gvcf
  Boolean output_vcf
  ## BQSR intervals
  String bqsr_intervals
  ## Readwriter read_filter args, starting after the QualCalFilter table
  String readwriter_readfilter_args
  ## Variant calling algorithm
  String calling_algo
  String calling_intervals
  ## Extra driver parameters
  String lc_driver_args
  String dedup_driver_args
  String readwriter_driver_args
  String bqsr_driver_args
  String calling_driver_args
  String genotyping_driver_args
  ## Extra algo parameters
  String bwa_args
  String bwa_chunk_size
  String sort_args
  String lc_args
  String dedup_args
  String readwriter_args
  String bqsr_args
  String calling_args
  String genotyping_args
  ## Alignment file formats
  Boolean alignment_cram
  Boolean dedup_cram
  Boolean readwriter_cram

  # Sentieon License configuration
  String sentieon_license_server
  Boolean use_instance_metadata
  String? sentieon_auth_mech
  String? sentieon_license_key

  # Execution configuration
  String disk
  String threads
  String memory
  Int preemptible_tries
  String docker
  String sentieon_release_dir
  
  # Some preprocessing
  Boolean make_gvcf = upload_gvcf || output_gvcf
  Boolean call_variants = make_gvcf || output_vcf
  Boolean run_genotyper
  String readwriter_suffix = if readwriter_cram then "cram" else "bam"
  String readwriter_index_suffix = if readwriter_cram then "cram.crai" else "bam.bai"
  String dollar = "$"
  command <<<
    set -exo pipefail
    mkdir -p /tmp
    export TMPDIR=/tmp

    # Check that the configuration is valid.
    # Supported variant callers are Genotyper, Haplotyper and DNAscope
    if [[ "${calling_algo}" != "Genotyper" && "${calling_algo}" != "Haplotyper" && "${calling_algo}" != "DNAscope" ]]; then
      echo "${calling_algo} is not a supported variant caller. Please set calling_algo to 'Genotyper', 'Haplotyper' or 'DNAscope'" >&2
      exit 1
    fi

    # Number of readgroups must match the number of fastq files
    first_fastq=(${sep=" " file_fastq_r1s})
    second_fastq=(${sep=" " file_fastq_r2s})
    read_groups=('${sep="' '" read_groups}')
    #Modified line, originally it was the same as the next one
    if [[ ${dollar}{#first_fastq[@]} -ne ${dollar}{#second_fastq[@]} ]]; then
      echo "The number of fastq files for r1 does not equal the number for r2"
      exit 1
    fi
    if [[ ${dollar}{#first_fastq[@]} -ne ${dollar}{#read_groups[@]} ]]; then
      echo "The number of fastq pairs does not equal the number of supplied readgroups"
      exit 1
    fi

    wait_list=()
    # License server setup
    license_file=""
    if [[ -n "$license_file" ]]; then
      # Using a license file
      echo "using license"
    elif [[ -n '${true="yes" false="" use_instance_metadata}' ]]; then
      python /opt/sentieon/gen_credentials.py ~/credentials.json ${default="''" sentieon_license_key} &
      sleep 5
      export SENTIEON_LICENSE=${default="" sentieon_license_server}
      export SENTIEON_AUTH_MECH=${default="" sentieon_auth_mech}
      export SENTIEON_AUTH_DATA=~/credentials.json
      read -r SENTIEON_JOB_TAG < ~/credentials.json.project
      export SENTIEON_JOB_TAG
    else
      export SENTIEON_LICENSE=${default="" sentieon_license_server}
      export SENTIEON_AUTH_MECH=${default="" sentieon_auth_mech}
    fi

    # Optimizations
    mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')
    export bwt_max_mem="$((mem_kb / 1024 / 1024 - 2))g"
    export MALLOC_CONF=lg_dirty_mult:-1

    # Alignment with BWA
    aligned_data=()
    for i in $(seq 1 ${dollar}{#first_fastq[@]}); do
      i=$((i - 1))
      LD_PRELOAD=${sentieon_release_dir}/lib/libjemalloc.so sentieon bwa mem -t ${threads} ${bwa_args} -K ${bwa_chunk_size} -R "${dollar}{read_groups[$i]}" ${ref_fasta} <(gsutil cp ${dollar}{first_fastq[$i]} -) <(gsutil cp ${dollar}{second_fastq[$i]} -) | \
        samblaster --addMateTags -a | \
        LD_PRELOAD=${sentieon_release_dir}/lib/libjemalloc.so sentieon util sort ${sort_args} -t ${threads} -i - --sam2bam -o ${sample_name}_sorted_${dollar}{i}.${true="cram" false="bam" alignment_cram}
      aligned_data+=(${sample_name}_sorted_${dollar}{i}.${true="cram" false="bam" alignment_cram})
    done

    export LD_PRELOAD=${sentieon_release_dir}/lib/libjemalloc.so

    # Dedup
    bam_input=""
    for f in ${dollar}{aligned_data[@]}; do
      bam_input="$bam_input -i $f"
    done
    sentieon driver ${lc_driver_args} -t ${threads} $bam_input -r ${ref_fasta} --algo LocusCollector ${lc_args} ${sample_name}_score.txt --algo MeanQualityByCycle ${sample_name}_mq_metrics.txt --algo QualDistribution ${sample_name}_qd_metrics.txt --algo GCBias --summary ${sample_name}_gc_summary.txt ${sample_name}_gc_metrics.txt --algo AlignmentStat --adapter_seq '' ${sample_name}_aln_metrics.txt --algo InsertSizeMetricAlgo ${sample_name}_is_metrics.txt
    sentieon driver ${dedup_driver_args} -t ${threads} $bam_input --algo Dedup ${dedup_args} --metrics ${sample_name}_dedup_metrics.txt --score_info ${sample_name}_score.txt --output_dup_read_name ${sample_name}_dup_qname.txt
    sentieon driver ${dedup_driver_args} -t ${threads} $bam_input --algo Dedup ${dedup_args} --dup_read_name ${sample_name}_dup_qname.txt ${sample_name}_deduped.${true="cram" false="bam" dedup_cram}
    for f in ${dollar}{aligned_data[@]}; do
      rm "$f" &
    done
    
    # Plot metrics - no need to add to $wait_list, runtime for each is only a few seconds
    sentieon plot GCBias -o ${sample_name}_gc-report.pdf ${sample_name}_gc_metrics.txt &
    sentieon plot QualDistribution -o ${sample_name}_qd-report.pdf ${sample_name}_qd_metrics.txt &
    sentieon plot MeanQualityByCycle -o ${sample_name}_mq-report.pdf ${sample_name}_mq_metrics.txt &
    sentieon plot InsertSizeMetricAlgo -o ${sample_name}_is-report.pdf ${sample_name}_is_metrics.txt &
    wait_list+=($!)


    sentieon driver ${"--interval " + bqsr_intervals} -r ${ref_fasta} -t ${threads} -i ${sample_name}_deduped.${true="cram" false="bam" dedup_cram} ${bqsr_driver_args} --algo QualCal -k ${sep=" -k " bqsr_vcfs} ${bqsr_args} ${sample_name}_recal.table
    
    # ReadWriter
    sentieon driver -r ${ref_fasta} -t ${threads} -i ${sample_name}_deduped.${true="cram" false="bam" dedup_cram} --read_filter QualCalFilter,table=${sample_name}_recal.table${readwriter_readfilter_args} ${readwriter_driver_args} --algo ReadWriter ${readwriter_args} ${sample_name}_recal.${readwriter_suffix}
    rm ${sample_name}_deduped.${true="cram" false="bam" dedup_cram} &

    # Ensure all output files are present so Cromwell does not error if they are streamed
    touch ${sample_name}_${calling_algo}.g.vcf.gz ${sample_name}_${calling_algo}.g.vcf.gz.tbi ${sample_name}_${calling_algo}.recal.vcf.gz ${sample_name}_${calling_algo}.recal.vcf.gz.tbi ${sample_name}.SNP.VQSR.pdf
    
    if [[ -n '${true="y" false="" call_variants}' ]]; then
      # Call variants
      sentieon driver -r ${ref_fasta} -t ${threads} -i ${sample_name}_recal.${readwriter_suffix} ${calling_driver_args} ${"--interval " + calling_intervals} --algo ${calling_algo} ${"-d " + dbsnp_vcf} --emit_conf 10 --call_conf 10 ${calling_args} ${true="--emit_mode GVCF" false="" make_gvcf} ${sample_name}_${calling_algo}${true=".g" false="" make_gvcf}.vcf.gz

      # this step is here to individually genotype samples - run_genotyper should be set to false if the aim is to joint genotype samples
      if [[ -n '${true="y" false="" run_genotyper}' ]]; then
        # Genotype the GVCF
        sentieon driver -r ${ref_fasta} -t ${threads} ${genotyping_driver_args} --algo GVCFtyper ${genotyping_args} --emit_conf 10 --call_conf 10 ${sample_name}_${calling_algo}.vcf.gz ${sample_name}_${calling_algo}.g.vcf.gz

        resource_text_SNP="--resource ${hapmap_vcf} --resource_param HapMap,known=false,training=true,truth=true,prior=15.0 "
        resource_text_SNP="$resource_text_SNP --resource ${omni_vcf} --resource_param Omni,known=false,training=true,truth=true,prior=12.0 "
        resource_text_SNP="$resource_text_SNP --resource ${onekg_vcf} --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
        resource_text_SNP="$resource_text_SNP --resource ${dbsnp_vcf} --resource_param dbSNP,known=true,training=false,truth=false,prior=2.0"

        resource_text_indel="--resource ${mills_vcf} --resource_param Mills,known=false,training=true,truth=true,prior=12.0 "
        resource_text_indel="$resource_text_indel --resource ${dbsnp_vcf} --resource_param dbSNP,known=true,training=false,truth=false,prior=2.0"
        resource_text_indel="$resource_text_indel --resource ${axiom_poly_vcf} --resource_param axiomPoly,known=false,training=true,truth=false,prior=10.0"

        
        #SNP RECAL
        sentieon driver -t 64 -r ${ref_fasta} --algo VarCal -v ${sample_name}_${calling_algo}.vcf.gz --max_gaussians 6 --tranches_file ${sample_name}.SNP.tranches  --var_type SNP --plot_file ${sample_name}.SNP.varcal.plot $resource_text_SNP --annotation QD --annotation MQ --annotation MQRankSum --annotation ReadPosRankSum --annotation FS --tranche 100 --tranche 99.95 --tranche 99.9 --tranche 99.8 --tranche 99.6 --tranche 99.5 --tranche 99.4 --tranche 99.3 --tranche 99.0 --tranche 98.0 --tranche 97.0 --tranche 90.0 ${sample_name}_${calling_algo}.vcf.SNP.recal

        sentieon driver -t 64 -r ${ref_fasta} --algo ApplyVarCal --sensitivity 99.7 -v ${sample_name}_${calling_algo}.vcf.gz --var_type SNP --recal ${sample_name}_${calling_algo}.vcf.SNP.recal --tranches_file ${sample_name}.SNP.tranches ${sample_name}_${calling_algo}.vcf.SNP.recaled.vcf.gz


        #INDEL RECAL
        sentieon driver -t 64 -r ${ref_fasta} --algo VarCal -v ${sample_name}_${calling_algo}.vcf.SNP.recaled.vcf.gz $resource_text_indel --max_gaussians 4 --annotation QD --annotation ReadPosRankSum --annotation FS --tranche 100.0 --tranche 99.95 --tranche 99.9 --tranche 99.5 --tranche 99.0 --tranche 97.0 --tranche 96.0 --tranche 95.0 --tranche 94.0 --tranche 93.5 --tranche 93.0 --tranche 92.0 --tranche 91.0 --tranche 90.0 --var_type INDEL --tranches_file ${sample_name}.INDEL.tranches ${sample_name}_${calling_algo}.vcf.SNP.INDEL.recal

        sentieon driver -t 64 -r ${ref_fasta} --algo ApplyVarCal -v ${sample_name}_${calling_algo}.vcf.SNP.recaled.vcf.gz --sensitivity 99.7 --var_type INDEL --recal ${sample_name}_${calling_algo}.vcf.SNP.INDEL.recal  --tranches_file ${sample_name}.INDEL.tranches ${sample_name}_${calling_algo}.recal.vcf.gz

        
       sentieon plot vqsr -o ${sample_name}.SNP.VQSR.pdf ${sample_name}.SNP.varcal.plot 
      fi    
    fi

    # Wait
    for pid in ${dollar}{wait_list[@]}; do
      wait $pid
    done

  >>>
  runtime {
    preemptible: preemptible_tries
    docker: docker
    memory: memory
    cpu: threads
    disks: disk
  }
  output {
    File alignment = "${sample_name}_recal.${readwriter_suffix}"
    File alignment_index = "${sample_name}_recal.${readwriter_index_suffix}"
    File gvcf = "${sample_name}_${calling_algo}.g.vcf.gz"
    File gvcf_index = "${sample_name}_${calling_algo}.g.vcf.gz.tbi"
    File dedup_metrics = "${sample_name}_dedup_metrics.txt"
    File mq_metrics = "${sample_name}_mq_metrics.txt"
    File mq_plot = "${sample_name}_mq-report.pdf"
    File qd_metrics = "${sample_name}_qd_metrics.txt"
    File qd_plot = "${sample_name}_qd-report.pdf"
    File gc_summary = "${sample_name}_gc_summary.txt"
    File gc_metrics = "${sample_name}_gc_metrics.txt"
    File gc_plot = "${sample_name}_gc-report.pdf"
    File as_metrics = "${sample_name}_aln_metrics.txt"
    File is_metrics = "${sample_name}_is_metrics.txt"
    File is_plot = "${sample_name}_is-report.pdf"
  }
}

# this step produces CCDG-compliant CRAMs
task stripTags {
    File cram
    File cram_index
    String sample_name
    
    File ref_fasta
    File ref_fasta_index
    
    String docker
    Int disk_size = ceil(size(cram, "GB")*4 + 50)
    
    command {
        mv ${cram} ./${sample_name}.cram
        mv ${cram_index} ./${sample_name}.cram.crai

        samtools view \
            ${sample_name}.cram \
            -h \
            -T ${ref_fasta} \
            -x XS \
            -C \
            > ${sample_name}_recal.cram
        
        samtools index ${sample_name}_recal.cram
    }
    
    output {
        File stripped_cram = "${sample_name}_recal.cram"
        File stripped_cram_index = "${sample_name}_recal.cram.crai"
    }
    
    runtime {
        docker: docker
        cpu: 1
        memory: "3.75 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
    }
}

task cramStats {
  File cram
  File cram_index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String sample_name

  String docker
  Int disk_size = ceil(size(cram, "GB")*2 + size(ref_fasta, "GB") + 10)

  command {
    samtools stats ${cram} --reference ${ref_fasta} > ${sample_name}_stats.txt
  }

  output {
    File stats = "${sample_name}_stats.txt"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 4
  }
}

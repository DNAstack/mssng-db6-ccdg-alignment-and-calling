workflow jointGenotype {
	File gvcf_URLs
	File region
	String joint_samplename

	File google_application_credentials_file

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
	File axiom_poly_vcf = "gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
	File axiom_poly_vcf_index = "gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"
 
	## Variant calling algorithm
	String calling_algo = "Haplotyper"
	## Extra driver parameters
	String genotyping_driver_args = ""
	## Extra algo parameters
	String genotyping_args = ""
	## Alignment file formats

	# Sentieon License configuration
	File? sentieon_license_file
	String sentieon_license_server = ""
	Boolean use_instance_metadata = false
	String? sentieon_auth_mech
	String? sentieon_license_key

	# Execution configuration
	String threads = "64"
	String memory = "240 GB"
	Int preemptible_tries = 1
	String sentieon_version = "201808.06"
	String docker = "dnastack/sentieon-bcftools:${sentieon_version}"


	call GVCFtyper {
		input:
			region = region,
			gvcf_URLs = gvcf_URLs,
			joint_samplename = joint_samplename,
			google_application_credentials_file = google_application_credentials_file,
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
			## Variant calling parameters
			calling_algo = calling_algo,
			## Extra driver parameters
			genotyping_driver_args = genotyping_driver_args,
			## Extra algo parameters
			genotyping_args = genotyping_args,
			# Sentieon License configuration
			sentieon_license_server = sentieon_license_server,
			sentieon_license_file = sentieon_license_file,
			use_instance_metadata = use_instance_metadata,
			sentieon_auth_mech = sentieon_auth_mech,
			sentieon_license_key = sentieon_license_key,
			# Execution configuration
			threads = threads,
			memory = memory,
			preemptible_tries = preemptible_tries,
			docker = docker
	}

	output {
		File partial_vcf = GVCFtyper.partial_vcf
		File partial_vcf_index = GVCFtyper.partial_vcf_index
	}
}

task GVCFtyper {
	File region
	String shard = basename(region, ".txt")
	File gvcf_URLs
	String joint_samplename

	File google_application_credentials_file

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

	## Variant calling algorithm
	String calling_algo
	## Extra driver parameters
	String genotyping_driver_args
	## Extra algo parameters
	String genotyping_args

	# Sentieon License configuration
	File? sentieon_license_file
	String sentieon_license_server
	Boolean use_instance_metadata
	String? sentieon_auth_mech
	String? sentieon_license_key

	# Execution configuration
	String threads
	String memory
	Int preemptible_tries
	String docker

	command <<<
		set -exo pipefail
		mkdir -p /tmp
		export TMPDIR=/tmp

		ulimit -s 327680

		# Check that the configuration is valid.
		# Supported variant callers are Genotyper, Haplotyper and DNAscope
		if [[ "${calling_algo}" != "Genotyper" && "${calling_algo}" != "Haplotyper" && "${calling_algo}" != "DNAscope" ]]; then
		  echo "${calling_algo} is not a supported variant caller. Please set calling_algo to 'Genotyper', 'Haplotyper' or 'DNAscope'" >&2
		  exit 1
		fi


		# License server setup
		license_file=${default="" sentieon_license_file}
		if [[ -n "$license_file" ]]; then
		  # Using a license file
		  export SENTIEON_LICENSE=${default="" sentieon_license_file}
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
		export VCFCACHE_BLOCKSIZE=4096
		export LD_PRELOAD=/opt/sentieon/sentieon-genomics-201808.06/lib/libjemalloc.so.1
		export MALLOC_CONF=lg_dirty_mult:-1

		# create .bed file for bcftools
		cat ${region} | tr " " "\n" | sed 's/\(.*\)-/\1\t/' | sed 's/\(.*\):/\1\t/' > ${shard}.bed

		export GOOGLE_APPLICATION_CREDENTIALS=${google_application_credentials_file}
        # Download the specified region from each gvcf file
        cat ${gvcf_URLs} | parallel --env GOOGLE_APPLICATION_CREDENTIALS \
                -j+0 \
                'export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token); \
                sample=$(basename {} _Haplotyper.g.vcf.gz); \
                bcftools view -R ${shard}.bed -Oz -o $sample\_${shard}.g.vcf.gz {}; \
                sentieon util vcfindex $sample\_${shard}.g.vcf.gz'

		# generate --shard string
		shard_command=""
		for region in $(cat ${region} | tr " " "\n"); do
			shard_command=$shard_command" --shard $region"
		done

		ls *.g.vcf.gz > input_files.txt
        
		# Genotype the GVCFs
		sentieon driver \
			-r ${ref_fasta} \
			-t ${threads} \
			${genotyping_driver_args} \
			$shard_command \
			--traverse_param 10000/200 \
			--algo GVCFtyper \
			--emit_conf 10 \
			--call_conf 10 \
      		--genotype_model multinomial \
			${genotyping_args} \
			${joint_samplename}_${shard}.vcf.gz \
			- < input_files.txt
	>>>

	output {
		File partial_vcf = "${joint_samplename}_${shard}.vcf.gz"
		File partial_vcf_index = "${joint_samplename}_${shard}.vcf.gz.tbi"
	}

	runtime {
		docker: docker
		cpu: threads
		memory: memory
		disks: "local-disk 4000 HDD"
		preemptible: preemptible_tries
	}
}

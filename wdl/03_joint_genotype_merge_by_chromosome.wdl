workflow mergeShards {
	Array [File] partial_vcfs
	Array [File] partial_vcf_indices
	String joint_samplename
	File gvcf_URLs
	# a .bed file containing the desired regions only
	File region

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

	# Sentieon License configuration
	File? sentieon_license_file
	String sentieon_license_server = ""
	Boolean use_instance_metadata = false
	String? sentieon_auth_mech
	String? sentieon_license_key

	# Execution configuration
	String threads = "4"
	String memory = "15 GB"
	String sentieon_version = "201808.06"
	String docker = "dnastack/sentieon-bcftools:${sentieon_version}"


	call mergePartialVCFs {
		input:
			partial_vcfs = partial_vcfs,
			partial_vcf_indices = partial_vcf_indices,
			joint_samplename = joint_samplename,
			gvcf_URLs = gvcf_URLs,
			region = region,
			# Known sites
			dbsnp_vcf = dbsnp_vcf,
			dbsnp_index = dbsnp_index,
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
			# Sentieon License configuration
			sentieon_license_server = sentieon_license_server,
			sentieon_license_file = sentieon_license_file,
			use_instance_metadata = use_instance_metadata,
			sentieon_auth_mech = sentieon_auth_mech,
			sentieon_license_key = sentieon_license_key,
			# Execution configuration
			threads = threads,
			memory = memory,
			docker = docker
	}

	output {
		File GVCFtyper_main_vcf = mergePartialVCFs.GVCFtyper_main_vcf
		File GVCFtyper_main_vcf_index = mergePartialVCFs.GVCFtyper_main_vcf_index
		File GVCFtyper_split_vcf = mergePartialVCFs.GVCFtyper_split_vcf
		File GVCFtyper_split_vcf_index = mergePartialVCFs.GVCFtyper_split_vcf_index
		File split_conf = mergePartialVCFs.split_conf
	}
}

task mergePartialVCFs {
	Array [File] partial_vcfs
	Array [File] partial_vcf_indices
	String joint_samplename
	File gvcf_URLs
	File region
	String chromosome = basename(region)

	# Known sites
	File? dbsnp_vcf
	File? dbsnp_index

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

	# Sentieon License configuration
	File? sentieon_license_file
	String sentieon_license_server
	Boolean use_instance_metadata
	String? sentieon_auth_mech
	String? sentieon_license_key

	# Execution configuration
	String threads
	String memory
	String docker


	command {
		set -exo pipefail
		mkdir -p /tmp
		export TMPDIR=/tmp

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
		export MALLOC_CONF=lg_dirty_mult:-1

		mv ${sep=' ' partial_vcfs} ${sep=' ' partial_vcf_indices} .

		# Generate split.conf file
		for line in $(cat ${gvcf_URLs})
		do
			sample=$(basename $line _Haplotyper.g.vcf.gz)
			echo $sample >> all_samples.txt
		done

		# we want all samples in a single file, and a 'main' file that can undergo VQSR
		split -d -l 100000 all_samples.txt group
		
		samples=$(cat group00 | sed 's/$/\t/g' | tr -d '\n' | sed 's/\t$//')
		echo -e "${joint_samplename}_GVCFtyper_file.vcf.gz\t$samples" >> split.conf

		# Merge GVCFs
		sentieon driver \
			--passthru \
			--algo GVCFtyper \
			--split_by_sample split.conf \
			--merge \
			${joint_samplename}_GVCFtyper_main.vcf.gz \
			$(ls *.vcf.gz)

		# Extract only the desired region from the main and samples files
		bcftools view \
			-R ${region} \
			-O z \
			-o ${joint_samplename}_GVCFtyper_file_${chromosome}.vcf.gz \
			${joint_samplename}_GVCFtyper_file.vcf.gz

		bcftools view \
			-R ${region} \
			-O z \
			-o ${joint_samplename}_GVCFtyper_main_${chromosome}.vcf.gz \
			${joint_samplename}_GVCFtyper_main.vcf.gz

		sentieon util vcfindex ${joint_samplename}_GVCFtyper_file_${chromosome}.vcf.gz
		sentieon util vcfindex ${joint_samplename}_GVCFtyper_main_${chromosome}.vcf.gz
	}

	output {
		File GVCFtyper_main_vcf = "${joint_samplename}_GVCFtyper_main_${chromosome}.vcf.gz"
		File GVCFtyper_main_vcf_index = "${joint_samplename}_GVCFtyper_main_${chromosome}.vcf.gz.tbi"
		File GVCFtyper_split_vcf = "${joint_samplename}_GVCFtyper_file_${chromosome}.vcf.gz"
		File GVCFtyper_split_vcf_index = "${joint_samplename}_GVCFtyper_file_${chromosome}.vcf.gz.tbi"
		File split_conf = "split.conf"
	}

	# No preemptible; will take > 24 hours
	runtime {
		docker: docker
		cpu: threads
		memory: memory
		disks: "local-disk 4000 HDD"
	}

}
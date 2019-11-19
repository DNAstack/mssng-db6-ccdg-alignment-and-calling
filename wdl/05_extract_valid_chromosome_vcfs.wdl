workflow extract {
	File vcf_main
	File vcf_main_index
	File vcf_samples
	File vcf_samples_index

	File all_samples

	String joint_samplename
	String region

	Boolean is_recalibrated

	# Sentieon License configuration
	File? sentieon_license_file
	String sentieon_license_server = ""
	Boolean use_instance_metadata = false
	String? sentieon_auth_mech
	String? sentieon_license_key

	# Execution configuration
	String threads = "2"
	String memory = "7.5 GB"
	String sentieon_version = "201808.06"
	String docker = "dnastack/sentieon-bcftools:${sentieon_version}"


	call extractSamples {
		input:
			vcf_main = vcf_main,
			vcf_main_index = vcf_main_index,
			vcf_samples = vcf_samples,
			vcf_samples_index = vcf_samples_index,
			all_samples = all_samples,
			joint_samplename = joint_samplename,
			region = region,
			is_recalibrated = is_recalibrated,
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
}

task extractSamples {
	File vcf_main
	File vcf_main_index
	File vcf_samples
	String vcf_samples_base = basename(vcf_samples)
	File vcf_samples_index

	File all_samples

	String joint_samplename
	String region

	# Only affects the name of the output
	Boolean is_recalibrated
	String recal = if is_recalibrated then ".recal" else ""

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
	Int disk_size = ceil((size(vcf_main, "GB") + size(vcf_samples, "GB"))*4 + 100)

	command {
		set -exo pipefail
		mkdir -p /tmp
		export TMPDIR=/tmp

		ulimit -s 327680

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

		mv ${vcf_samples} ${vcf_samples_index} .

		# Regenerate split.conf
		echo -ne "${vcf_samples_base}\t" | cat - ${all_samples} > split.conf

		# Generate samples.csv
		cat ${all_samples} | tr '\t' , > samples.csv

		extract.sh \
			${vcf_main} \
			split.conf \
			$(cat samples.csv) \
			| bgzip -@ ${threads} > ${joint_samplename}.${region}${recal}.vcf.gz

		sentieon util vcfindex ${joint_samplename}.${region}${recal}.vcf.gz
	}

	output {
		File vcf = "${joint_samplename}.${region}${recal}.vcf.gz"
		File vcf_index = "${joint_samplename}.${region}${recal}.vcf.gz.tbi"
	}

	# no preemptible; will run >24h
	runtime {
		docker: docker
		cpu: threads
		memory: memory
		disks: "local-disk " + disk_size + " HDD"
	}
}
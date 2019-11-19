workflow variantMetrics {
	File vcf
	File vcf_index
	String joint_sample_name
	String region

	File ref_dict = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict"
	File dbsnp_vcf = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
	File dbsnp_vcf_index = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

	String docker = "dnastack/picard_samtools:2.18.9"

	call runVariantMetrics {
		input:
			vcf = vcf,
			vcf_index = vcf_index,
			joint_sample_name = joint_sample_name,
			region = region,
			ref_dict = ref_dict,
			dbsnp_vcf = dbsnp_vcf,
			dbsnp_vcf_index = dbsnp_vcf_index,
			docker = docker
	}

	output {
		File summary_metrics_file = runVariantMetrics.summary_metrics_file
		File detail_metrics_file = runVariantMetrics.detail_metrics_file
	}

}

task runVariantMetrics {
	File vcf
	File vcf_index
	String joint_sample_name
	String region

	File ref_dict
	File dbsnp_vcf
	File dbsnp_vcf_index

	String docker
	Int disk_size = ceil((size(vcf, "GB") + size(dbsnp_vcf, "GB")) * 2 + 50)

	command {
	    java -Xmx6g -Xms6g -jar $PICARD \
			CollectVariantCallingMetrics \
			INPUT=${vcf} \
			DBSNP=${dbsnp_vcf} \
			SEQUENCE_DICTIONARY=${ref_dict} \
			OUTPUT=${joint_sample_name}.${region} \
			THREAD_COUNT=2
	}

	output {
    	File summary_metrics_file = "${joint_sample_name}.${region}.variant_calling_summary_metrics"
    	File detail_metrics_file = "${joint_sample_name}.${region}.variant_calling_detail_metrics"
  	}

	runtime {
		docker: docker
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk " + disk_size + " HDD"
	}
}

version 1.0

workflow msisensor {
	input {
		File normalbam
		File normalbai
		File tumorbam
		File tumorbai
		String basename = basename("~{tumorbam}", ".bam")
		Boolean boostrapIt
	}

	parameter_meta {
		normalbam: "normal input .bam file"
		normalbai: "normal input .bai file"
		tumorbam: "tumor input .bam file"
		tumorbai: "tumor input .bai file"
		basename: "Base name"
	}

	call runMSIsensor {
		input: 
			normalbam = normalbam,
			tumorbam = tumorbam,
			normalbai = normalbai,
			tumorbai = tumorbai
	}

	if(boostrapIt == true){

		call bootstrapMSIsensor {
			input:
				normalbam = normalbam,
				tumorbam = tumorbam,
				normalbai = normalbai,
				tumorbai = tumorbai
		}

	}

	meta {
		author: "Felix Beaudry"
		email: "fbeaudry@oicr.on.ca"
		description: "Microsatelite Instability (MSI) detection using msisensor-pro"
		dependencies: 
		[
			{
				name: "msisensorpro/1.2.0",
				url: "https://github.com/broadinstitute/gatk/releases"
			}
		]
		output_meta: {
			msiFinalOutput: "Final msisensor call as .tsv, last column is msi score",
			msibooted: "msisensor calls bootstrapped"
		}
	}
	output {
		File msiFinalOutput = "~{basename}.msi"
		File msibooted = "~{basename}.msi.booted"
	}
}

task runMSIsensor {
	input {
		File normalbam 
		File tumorbam
		File normalbai 
		File tumorbai 
		String basename = basename("~{tumorbam}", ".bam")
		String modules = "msisensorpro/1.2.0 msisensor-microsatlist/hg38p12"
		String msifile = "$MSISENSOR_MICROSATLIST_ROOT/hg38_random.fa.list"
		String? difficultRegions
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		normalbam: "normal input .bam file"
		tumorbam: "tumor input .bam file"
		normalbai: "normal input .bai file"
		tumorbai: "tumor input .bai file"
		basename: "Base name"
		modules: "Required environment modules"
		msifile: "list of microsats identified by msisensor-scan"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
		difficultRegions: "Path to .bed of difficult regions to exclude, string must include the flag -e "
	}

	command <<<
		set -euo pipefail

		msisensor-pro msi \
			-d ~{msifile} ~{difficultRegions} \
			-n ~{normalbam} -t ~{tumorbam} \
			-o ~{basename}.msi 

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File msiGermline = "~{basename}.msi_germline"
		File msiSomatic = "~{basename}.msi_somatic"
		File msiDistribution = "~{basename}.msi_dis"
		File? msiFinalOutput = "~{basename}.msi"
	}

	meta {
		output_meta: {
			msiGermline: "MSI calls for germline",
			msiSomatic: "MSI calls for soma",
			msiDistribution: "MSI distribution per site for normal and tumor",
			msiFinalOutput: "Final msisensor call as .tsv, last column is msi score"
		}
	}
}

task bootstrapMSIsensor {
	input {
		Int boots = 100
		Int loci = 500
		File normalbam 
		File tumorbam
		File normalbai 
		File tumorbai 
		String basename = basename("~{tumorbam}", ".bam")
		String modules = "msisensorpro/1.2.0 msisensor-microsatlist/hg38p12"
		String msifile = "$MSISENSOR_MICROSATLIST_ROOT/hg38_random.fa.list"
		String? difficultRegions
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
		String? difficultRegions
	}

	parameter_meta {
		boots: "number of bootstraps"
		loci: "number of loci to include in each bootstrap"
		normalbam: "normal input .bam file"
		tumorbam: "tumor input .bam file"
		normalbai: "normal input .bai file"
		tumorbai: "tumor input .bai file"
		basename: "Base name"
		modules: "Required environment modules"
		msifile: "list of microsats identified by msisensor-scan"
		difficultRegions: "bed file of regions to avoid, if necessary"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
		difficultRegions: "Path to .bed of difficult regions to exclude, string must include the flag -e "
	}

	command <<<
		set -euo pipefail

		for boot in {1..~{boots}}
		do
			shuf -n ~{loci} ~{msifile} >rep.list

			sort -k1,1 -k2,2n rep.list >rep.list.sorted

			msisensor-pro msi \
				-d rep.list.sorted ~{difficultRegions} \
				-n ~{normalbam} -t ~{tumorbam} \
				-o ~{basename}.msi

			awk -v boot="${boot}" '$1 !~ "Total_Number_of_Sites" {print boot"\t"$1"\t"$2"\t"$3}' ~{basename}.msi >>~{basename}.msi.booted

		done

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File msibooted = "~{basename}.msi.booted"
	}

	meta {
		output_meta: {
			msibooted: "MSI calls for germline, bootstrapped"
		}
	}
}


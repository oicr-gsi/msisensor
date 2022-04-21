version 1.0

workflow msisensor {
	input {
		File normalbam
		File tumorbam
		String basename = basename("~{tumorbam}", ".bam")
	}

	parameter_meta {
		normalbam: "normal input .bam file"
		tumorbam: "tumor input .bam file"
	}

	call msisensor {
		input: 
			normalbam = normalbam,
			tumorbam = tumorbam
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
			msiFinalOutput: "Final msisensor call as .tsv, last column is msi score"
		}
	}
	output {
		File msiFinalOutput = "~{basename}.msi"
	}
}

task msisensor {
	input {
		File normalbam 
		File tumorbam 
		String basename = basename("~{tumorbam}", ".bam")
		String modules = "msisensorpro/1.2.0 msisensor-microsatlist/hg38p12"
		String msifile = "${MSISENSOR_MICROSATLIST_ROOT}/hg38_random.fa.list"
		String? difficultRegions
		Int jobMemory = 5
		Int threads = 10
		Int timeout = 10
	}

	parameter_meta {
		normalbam: "normal input .bam file"
		tumorbam: "tumor input .bam file"
		basename: "Base name"
		modules: "Required environment modules"
		msifile: "list of microsats identified by msisensor-scan"
		difficultRegions: "bed file of regions to avoid, if necessary"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		msisensor-pro msi \
			-d ~{msifile} \
			-n ~{normalbam} -t ~{tumorbam} \
			-o ~{basename}.msi \
			-b ~{threads} 

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
		File msiFinalOutput = "~{basename}.msi"
	}

	meta {
		output_meta: {
			msiGermline: "MSI calls for germline"
			msiSomatic: "MSI calls for soma"
			msiDistribution: "MSI distribution per site for normal and tumor",
			msiFinalOutput: "Final msisensor call as .tsv, last column is msi score"
		}
	}
}

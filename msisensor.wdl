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

	call runMSIsensor as WG_MSIsensor {
		input: 
			normalbam = normalbam,
			tumorbam = tumorbam,
			normalbai = normalbai,
			tumorbai = tumorbai
	}

	if(boostrapIt == true){
		call make_boot_interval{}

		scatter ( boot in make_boot_interval.boot_interval ) {

			call make_boots {}

			call runMSIsensor as boot_MSIsensor {
				input: 
					normalbam = normalbam,
					tumorbam = tumorbam,
					normalbai = normalbai,
					tumorbai = tumorbai,
					msifile = make_boots.msilist_boot
			}

		}

		call gather_boots {
			input:
				msiFinalOutputs = select_all(boot_MSIsensor.msiFinalOutput)
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
			msiFinalOutput: "Final msisensor call as .tsv, last column is msi score"
		}
	}
	output {
		File msiFinalOutput = "~{basename}.msi"
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

task make_boots {
	input {
		Int loci = 5000
		String modules = "msisensorpro/1.2.0 msisensor-microsatlist/hg38p12"
		String msifile = "$MSISENSOR_MICROSATLIST_ROOT/hg38_random.fa.list"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		loci: "number of loci to include in each bootstrap"
		modules: "Required environment modules"
		msifile: "list of microsats identified by msisensor-scan"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		shuf -n ~{loci} ~{msifile} >rep.list

		sort -k1,1 -k2,2n rep.list >rep.list.sorted


	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File msilist_boot = "rep.list.sorted"
	}

	meta {
		output_meta: {
			msilist_boot: "random subset of MSI sites"
		}
	}
}

task make_boot_interval {
	input {
		String bootstraps = 1000
		Int jobMemory = 4
		Int timeout = 12
	}

	parameter_meta {
		bootstraps: "number of bootstraps"
		jobMemory: "Memory for this task in GB"
		timeout: "Timeout in hours, needed to override imposed limits"
	}

	command <<<
		python <<CODE
		import os, re

		for boot in range(~{bootstraps}):
			print(boot)
		CODE
	>>>

	runtime {
		memory:  "~{jobMemory} GB"
		timeout: "~{timeout}"
	}

	output {
		Array[Array[String]] boot_interval = read_tsv(stdout())
	}
}

task gather_boots {
	input {
		Array[File] msiFinalOutputs
		Int jobMemory = 4
		Int timeout = 12
	}

	parameter_meta {
		msiFinalOutputs: "msiFinalOutput files from bootstraps"
		jobMemory: "Memory for this task in GB"
		timeout: "Timeout in hours, needed to override imposed limits"
	}

	command <<<
		set -euo pipefail

		cat ~{sep=' ' msiFinalOutputs} >msiFinalOutput.booted.txt


	>>>

	runtime {
		memory:  "~{jobMemory} GB"
		timeout: "~{timeout}"
	}

	output {
		File booted_msi = "msiFinalOutput.booted.txt"
	}
}


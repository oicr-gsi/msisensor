# msisensor

Microsatelite Instability (MSI) detection using msisensor-pro

## Overview

## Dependencies

* [msisensorpro 1.2.0](https://github.com/broadinstitute/gatk/releases)


## Usage

### Cromwell
```
java -jar cromwell.jar run msisensor.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`normalbam`|File|normal input .bam file
`normalbai`|File|normal input .bai file
`tumorbam`|File|tumor input .bam file
`tumorbai`|File|tumor input .bai file


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`basename`|String|basename("~{tumorbam}",".bam")|Base name


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`runMSIsensor.basename`|String|basename("~{tumorbam}",".bam")|Base name
`runMSIsensor.modules`|String|"msisensorpro/1.2.0 msisensor-microsatlist/hg38p12"|Required environment modules
`runMSIsensor.msifile`|String|"$MSISENSOR_MICROSATLIST_ROOT/hg38_random.fa.list"|list of microsats identified by msisensor-scan
`runMSIsensor.difficultRegions`|String?|None|bed file of regions to avoid, if necessary
`runMSIsensor.jobMemory`|Int|64|Memory allocated for this job (GB)
`runMSIsensor.threads`|Int|1|Requested CPU threads
`runMSIsensor.timeout`|Int|10|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`msiFinalOutput`|File|Final msisensor call as .tsv, last column is msi score


## Commands
 This section lists the core command run by the msisensor workflow
 
 * Running runMSIsensor
 
 MSIsensor (in classic mode) requires paired tumor-normal .bam files (and .bai files) to produce an MSI score (fraction of MS sites showing instability, for which MSI is positive at scores higher than 3%. The command will also output the somatic/germline calls and the distribution of sites. 
 
 		msisensor-pro msi \
 			-d ~{msifile} \
 			-n ~{normalbam} -t ~{tumorbam} \
 			-o ~{basename}.msi \
 			-b ~{threads}  


## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_

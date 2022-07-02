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
`boostrapIt`|Boolean|Enable bootstrapping


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|basename("~{tumorbam}",".bam")|Base name


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`runMSIsensor.outputFileNamePrefix`|String|basename("~{tumorbam}",".bam")|Base name
`runMSIsensor.modules`|String|"msisensorpro/1.2.0 msisensor-microsatlist/hg38p12"|Required environment modules
`runMSIsensor.msifile`|String|"$MSISENSOR_MICROSATLIST_ROOT/hg38_random.fa.list"|list of microsats identified by msisensor-scan
`runMSIsensor.difficultRegions`|String?|None|Path to .bed of difficult regions to exclude
`runMSIsensor.jobMemory`|Int|64|Memory allocated for this job (GB)
`runMSIsensor.threads`|Int|4|Requested CPU threads
`runMSIsensor.timeout`|Int|10|Hours before task timeout
`bootstrapMSIsensor.boots`|Int|100|number of bootstraps
`bootstrapMSIsensor.loci`|Int|500|number of loci to include in each bootstrap
`bootstrapMSIsensor.outputFileNamePrefix`|String|basename("~{tumorbam}",".bam")|Base name
`bootstrapMSIsensor.modules`|String|"msisensorpro/1.2.0 msisensor-microsatlist/hg38p12"|Required environment modules
`bootstrapMSIsensor.msifile`|String|"$MSISENSOR_MICROSATLIST_ROOT/hg38_random.fa.list"|list of microsats identified by msisensor-scan
`bootstrapMSIsensor.jobMemory`|Int|64|Memory allocated for this job (GB)
`bootstrapMSIsensor.threads`|Int|4|Requested CPU threads
`bootstrapMSIsensor.timeout`|Int|10|Hours before task timeout
`bootstrapMSIsensor.difficultRegions`|String?|None|bed file of regions to avoid, if necessary


### Outputs

Output | Type | Description
---|---|---
`msiFinalOutput`|File|Final msisensor call as .tsv, last column is msi score
`msibooted`|File|msisensor calls bootstrapped


## Commands
 This section lists the core command run by the msisensor workflow
 
 * Running msisensor
 
 Microsatelite Instability (MSI) detection using msisensor-pro. MSIsensor-pro is an updated version of msisensor.
 MSIsensor-pro evaluates Microsatellite Instability (MSI) for cancer patients with NGS data.
 It accepts the whole genome sequencing, whole exome sequencing and target region (panel) sequencing data as input.
 MSIsensor-pro introduces a multinomial distribution model to quantify polymerase slippages for each tumor sample
 and a discriminative sites selection method to enable MSI detection without matched normal samples.
 
 Run msisensor-pro (the main step)
 
 '''
  set -euo pipefail
 
  msisensor-pro msi
    -d MSIFILE 
    -e difficultRegions (Optional)
    -n NORMAL_BAM 
    -t TUMOR_BAM
    -o BASENAME.msi 
 
 '''
 
 Optional bootstrapping with the same inputs:
 
 '''
  set -euo pipefail
 
  for boot in {1..~{boots}}
  do
  shuf -n ~{loci} ~{msifile} > rep.list
 
  sort -k1,1 -k2,2n rep.list > rep.list.sorted
 
  msisensor-pro msi 
   -d rep.list.sorted 
   -e difficultRegions (Optional)
   -n NORMAL_BAM
   -t TUMOR_BAM
   -o BASENAME.msi
 
  awk -v boot="${boot}" '$1 !~ "Total_Number_of_Sites" {print boot"\t"$1"\t"$2"\t"$3}' ~{basename}.msi >> ~{basename}.msi.booted
  done
 
 '''
 
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_

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
`reference`|String|reference genome of input sample


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|basename("~{tumorbam}",".bam")|Base name


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`runMSIsensor.outputFileNamePrefix`|String|basename("~{tumorbam}",".bam")|Base name
`runMSIsensor.difficultRegions`|String?|None|Path to .bed of difficult regions to exclude
`runMSIsensor.jobMemory`|Int|64|Memory allocated for this job (GB)
`runMSIsensor.threads`|Int|4|Requested CPU threads
`runMSIsensor.timeout`|Int|10|Hours before task timeout
`bootstrapMSIsensor.boots`|Int|100|number of bootstraps
`bootstrapMSIsensor.loci`|Int|500|number of loci to include in each bootstrap
`bootstrapMSIsensor.outputFileNamePrefix`|String|basename("~{tumorbam}",".bam")|Base name
`bootstrapMSIsensor.jobMemory`|Int|64|Memory allocated for this job (GB)
`bootstrapMSIsensor.threads`|Int|4|Requested CPU threads
`bootstrapMSIsensor.timeout`|Int|10|Hours before task timeout
`bootstrapMSIsensor.difficultRegions`|String?|None|bed file of regions to avoid, if necessary


### Outputs

Output | Type | Description | Labels
---|---|---|---
`msiGermline`|File|A poorly documented output, ostensibly germline-specific metrics for MSI sites|vidarr_label: msiGermline
`msiSomatic`|File|A poorly documented output, ostensibly somatic-specific metrics for MSI sites|vidarr_label: msiSomatic
`msiFinalOutput`|File|Final msisensor call as .tsv, last column is msi score|vidarr_label: msiFinalOutput
`msibooted`|File?|msisensor calls bootstrapped|vidarr_label: msibooted


## Commands
 
This section lists command(s) run by msisensor workflow
 
* Running msisensor
 
Microsatelite Instability (MSI) detection using msisensor-pro. MSIsensor-pro is an updated version of msisensor.
MSIsensor-pro evaluates Microsatellite Instability (MSI) for cancer patients with NGS data.
It accepts the whole genome sequencing, whole exome sequencing and target region (panel) sequencing data as input.
MSIsensor-pro introduces a multinomial distribution model to quantify polymerase slippages for each tumor sample
and a discriminative sites selection method to enable MSI detection without matched normal samples.
 
### Run msisensor-pro (the main step)
 
```
  set -euo pipefail
 
  msisensor-pro msi
    -d MSIFILE 
    -e difficultRegions (Optional)
    -n NORMAL_BAM 
    -t TUMOR_BAM
    -o BASENAME.msi 
 
```
 
### Optional bootstrapping with the same inputs:
 
```
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
 
```
 
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_

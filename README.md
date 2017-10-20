# TCGA_germline-nf
extraction of germline variants from TCGA data

#### Dependencies

Install [nextflow](http://www.nextflow.io/).

```bash
curl -fsSL get.nextflow.io | bash
```
And move it to a location in your `$PATH` (`/usr/local/bin` for example here):
```bash
sudo mv nextflow /usr/local/bin
```

#### Description

This program takes in input a TCGA data-like folder of VCF files, _e.g._ `myfolder` containing all `/ids_sample1/sample1.vcf`, `/ids_sample2/sample2.vcf`, `/ids_sample3/sample3.vcf`, etc.  
Then it filters out germline variants from these VCF files, with a minimum coverage `--min_DP` and variant allelic fraction `--min_af`.  
After that it reformats these files into annovar inputs (by cancer type) and output them in a result folder (`--out_folder`).  

#### Execution
Nextflow seamlessly integrates with GitHub hosted code repositories:

`nextflow run iarcbioinfo/TCGA_germline-nf --TCGA_folder myfolder`

#### Help and options
You can print the help manual by providing `--help` in the execution command line:
```bash
nextflow run iarcbioinfo/TCGA_germline-nf --help
```

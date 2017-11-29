# TCGA_platypus-nf
convert TCGA VCF into annovar-ready tables with nextflow

#### Dependencies

Install [nextflow](http://www.nextflow.io/).

```bash
curl -fsSL get.nextflow.io | bash
```
And move it to a location in your `$PATH` (`/usr/local/bin` for example here):
```bash
sudo mv nextflow /usr/local/bin
```  

Install [vt](https://genome.sph.umich.edu/wiki/Vt#Installation), and put the executable in your PATH.  

Install [VCFtools](https://github.com/vcftools/vcftools.), and put `vcf-sort` in your PATH.

#### Description

This program takes in input a folder of compressed VCF files, from variant calling on TCGA data : 10th column of header must be the TCGA sample barcode (see [specification here](https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode)).  
It reformats these VCF into annovar input tables by cancer type and output them in a result folder (`--out_folder`).
__Caution__: VCF must be zipped, with extension `vcf.gz`.  

#### Execution
Nextflow seamlessly integrates with GitHub hosted code repositories:

`nextflow run iarcbioinfo/TCGA_platypus-nf --TCGA_folder myfolder`

#### Help and options
You can print the help manual by providing `--help` in the execution command line:
```bash
nextflow run iarcbioinfo/TCGA_platypus-nf --help
```

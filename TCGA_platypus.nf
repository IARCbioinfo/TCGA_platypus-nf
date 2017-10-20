params.help = null
params.out_folder = "."

if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------------'
    log.info '                 CONVERT TCGA VCF TO ANNOVAR-READY TABLES            '
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run TCGA_platypus.nf --TCGA_folder /path/to/parent/TCGA/folder'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --TCGA_folder        FOLDER                  Input folder of TCGA VCFs.'
    log.info 'Optional arguments:'
    log.info '    --out_folder         FOLDER                  Output folder containing annovar-ready tables.'
    log.info ''
    log.info ''
    exit 1
}

all_vcf = Channel.fromPath( params.TCGA_folder+'/*.vcf.gz')
                 .ifEmpty { error "empty TCGA folders" }

process reformat {

  tag { SM_tag }

  input:
  file all_vcf

  output:
  file '*reformat.tsv' into reformated

  shell:
  SM_tag = all_vcf.baseName.substring(0,12)
  '''
  zcat !{all_vcf} > uncompressed.vcf
  reformat.sh uncompressed.vcf !{baseDir}/data/tissueSourceSite.tsv !{baseDir}/data/diseaseStudy.tsv
  rm uncompressed.vcf
  '''

}

process merge {

  publishDir params.out_folder, mode: 'move'

  input:
  file all_reformated from reformated.toList()

  output:
  file "*.tsv" into reformated_for_annovar mode flatten

  shell:
  '''
  cat *reformat.tsv > big.tsv
  awk -F" " '{print >  "TCGA_platypus_reformat_"$NF".tsv"}' big.tsv
  rm big.tsv
  '''
}

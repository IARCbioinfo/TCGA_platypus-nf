params.help = null
params.out_folder = "."
params.min_af = 0.1
params.min_DP = 10

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
    log.info '    --min_af             VALUE                   Minimum allelic fraction to consider a germline. Default=0.1.'
    log.info '    --min_DP             VALUE                   Minimum coverage to consider a site. Default=10.'
    log.info ''
    log.info ''
    exit 1
}

all_vcf = Channel.fromPath( params.TCGA_folder+'/*.vcf.gz')
                 .ifEmpty { error "empty TCGA folders" }

process germline_filter {

  tag { SM_tag }

  input:
  file all_vcf

  output:
  set val(SM_tag), file("*filter.vcf") into germ_filt

  shell:
  SM_tag = all_vcf.baseName.substring(0,12)
  '''
  zcat !{all_vcf} > uncompressed.vcf
  filter_germline.r --vcf=uncompressed.vcf --min_af=!{params.min_af} --min_DP=!{params.min_DP}
  '''

}

process reformat {

  tag { SM_tag }

  input:
  set val(SM_tag), file("*filter.vcf") from germ_filt

  output:
  file '*reformat.tsv' into reformated

  shell:
  '''
  reformat.sh filter.vcf !{baseDir}/data/tissueSourceSite.tsv !{baseDir}/data/diseaseStudy.tsv
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

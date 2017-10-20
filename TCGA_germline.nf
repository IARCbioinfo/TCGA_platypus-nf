params.help = null
params.min_af = 0.25
params.min_DP = 20
params.out_folder = "."

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '              TCGA GERMLINE EXTRACTION            '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run TCGA_germline.nf --min_af 0.25 --min_DP 20 --TCGA_folder /path/to/parent/TCGA/folder'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --TCGA_folder        FOLDER                  Parent folder of TCGA data.'
    log.info 'Optional arguments:'
    log.info '    --min_af             VALUE                   Minimum allelic fraction to consider a germline. Default=0.25.'
    log.info '    --min_DP             VALUE                   Minimumk coverage to consider a site. Default=20.'
    log.info ''
    log.info ''
    exit 1
}

all_vcf = Channel.fromPath( params.TCGA_folder+'/*/*.vcf')
                 .ifEmpty { error "empty TCGA subfolders" }

process germline_filter {

  tag { SM_tag }

  input:
  file all_vcf

  output:
  set val(SM_tag), file("*filter.vcf") into germ_filt

  shell:
  SM_tag = all_vcf.baseName.substring(0,12)
  '''
  filter_germline.r --vcf=!{all_vcf} --min_af=!{params.min_af} --min_DP=!{params.min_DP}
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
  awk -F" " '{print >  "TCGA_germline_reformat_"$NF".tsv"}' big.tsv
  rm big.tsv
  '''
}

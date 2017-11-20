params.help = null
params.out_folder = "."
params.min_af = 0.1
params.min_DP = 10
params.blood_tissue_filter = false
params.min_blood_QUAL = 100
params.min_tissue_QUAL = 100
params.cpu_R = 1
params.mem_R = 4


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
    log.info '    --TCGA_folder          FOLDER                  Input folder of TCGA VCFs.'
    log.info 'Optional arguments:'
    log.info '    --out_folder           FOLDER                  Output folder containing annovar-ready tables.'
    log.info '    --min_af               VALUE                   Minimum allelic fraction to consider a germline. Default=0.1.'
    log.info '    --min_DP               VALUE                   Minimum coverage to consider a site. Default=10.'
    log.info '    --min_blood_QUAL       VALUE                   Minimum QUAL to retain a variant found only in the blood. Default=100.'
    log.info '    --min_tissue_QUAL      VALUE                   Minimum QUAL to retain a variant found only in the tissue. Default=100.'
    log.info '    --cpu_R                VALUE                   Number of cpu to run R filtering on blood/normal samples in parallel. Default=1.'
    log.info '    --mem_R                VALUE                   Memory (in GB) for R filtering on blood/normal samples. Default=4.'
    log.info "Flags:"
    log.info '    --blood_tissue_filter                          To filter callings if both blood and tissue samples are available.'
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

process filter_blood_tissue {

  cpus params.cpu_R
  memory params.mem_R+'G'

  input:
  file all_reformated from reformated.toList()

  output:
  file "*blood_tissue_filtered.tsv" into blood_tissue_filtered

  shell:
  if(params.blood_tissue_filter){
    '''
    blood_tissue_filter.R --min_blood_QUAL=!{params.min_blood_QUAL} --min_tissue_QUAL=!{params.min_tissue_QUAL}
    '''
  } else {
    '''
      for file in *.tsv
      do
        ln -s "$file" "${file/.tsv/_blood_tissue_filtered.tsv}"
      done
    '''
  }
}

process merge {

  publishDir params.out_folder, mode: 'move'

  input:
  file all_filtered from blood_tissue_filtered

  output:
  file "*.tsv" into reformated_for_annovar mode flatten

  shell:
  '''
  cat *blood_tissue_filtered.tsv > big.tsv
  awk -F" " '{print >  "TCGA_platypus_reformat_"$NF".tsv"}' big.tsv
  rm big.tsv
  '''
}

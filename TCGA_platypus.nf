params.help = null
params.out_folder = "."
params.min_af = 0.1
params.min_DP = 10
params.blood_tissue_filter = false

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
    log.info '    --ref                FILE (with index)       Reference fasta file indexed.'
    log.info "Flags:"
    log.info '    --blood_tissue_filter                          To filter callings if both blood and tissue samples are available.'
    log.info ''
    log.info ''
    exit 1
}

assert (params.ref != true) && (params.ref != null) : "please specify --ref option (--ref reference.fasta(.gz))"
assert (params.TCGA_folder != true) && (params.TCGA_folder != null) : "please specify --TCGA_folder option"

fasta_ref = file(params.ref)
fasta_ref_fai = file( params.ref+'.fai' )

vcf = Channel.fromPath( params.TCGA_folder+'/*.vcf.gz')
                 .ifEmpty { error "empty TCGA folders" }

process vt {

  tag { vcf_tag }

  input:
  file vcf
  file fasta_ref
  file fasta_ref_fai

  output:
  file("${vcf_tag}_vt.vcf.gz") into vt_VCF

  shell:
  vcf_tag = vcf.baseName.replace(".vcf","")
  '''
  vcf-sort !{vcf_tag}.vcf.gz | vt decompose -s - | vt decompose_blocksub -a - | vt normalize -r !{fasta_ref} -q - | vt uniq - | bgzip -c > !{vcf_tag}_vt.vcf.gz
  '''
}

process germline_filter {

  tag { SM_tag }

  input:
  file vt_VCF

  output:
  set val(SM_tag), file("*filter.vcf") into germ_filt

  shell:
  SM_tag = vt_VCF.baseName.substring(0,12)
  '''
  zcat !{vt_VCF} > uncompressed.vcf
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
    blood_tissue_filter.R
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

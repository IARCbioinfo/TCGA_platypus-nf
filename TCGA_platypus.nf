params.help = null
params.out_folder = "."
params.min_af = 0.1
params.min_DP = 10
params.blood_tissue_filter = false
params.annovar_cpu = 1
params.annovar_db = null

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
    log.info '    --annovar_db           FOLDER                  Path to annovar database.'
    log.info 'Optional arguments:'
    log.info '    --out_folder           FOLDER                  Output folder containing annovar-ready tables.'
    log.info '    --min_af               VALUE                   Minimum allelic fraction to consider a germline. Default=0.1.'
    log.info '    --min_DP               VALUE                   Minimum coverage to consider a site. Default=10.'
    log.info '    --annovar_cpu          VALUE                   Number of threads used by annovar. Default=1'
    log.info "Flags:"
    log.info '    --blood_tissue_filter                          To filter callings if both blood and tissue samples are available.'
    log.info ''
    log.info ''
    exit 1
}

assert (params.TCGA_folder != true) && (params.TCGA_folder != null) : "please specify --TCGA_folder option"
assert (params.annovar_db != true) && (params.annovar_db != null) : "please specify --annovar_db option"

vcf = Channel.fromPath( params.TCGA_folder+'/*.vcf.gz')
                 .ifEmpty { error "empty TCGA folders" }

process convert2annovar {

  tag { SM_tag }

  input:
  file vcf

  output:
  file("${vcf_tag}_convert.vcf") into convert_VCF

  shell:
  SM_tag = vcf.baseName.substring(0,12)
  vcf_tag = vcf.baseName.replace(".vcf","")
  '''
  zcat !{vcf} | sed '/^#CHROM/Q' > !{vcf_tag}_convert.vcf
  zcat !{vcf} | grep -m1 "^#CHROM" | sed 's/POS\tID/START\tEND/g' >> !{vcf_tag}_convert.vcf
  convert2annovar.pl -format vcf4 -includeinfo !{vcf} |  grep -v "^#" | cut -f-5,11- >> !{vcf_tag}_convert.vcf
  '''

}

process germline_filter {

  tag { SM_tag }

  input:
  file convert_VCF

  output:
  set val(SM_tag), file("*filter.vcf") into germ_filt

  shell:
  SM_tag = convert_VCF.baseName.substring(0,12)
  '''
  filter_germline.r --vcf=!{convert_VCF} --min_af=!{params.min_af} --min_DP=!{params.min_DP}
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

  input:
  file all_reformated from reformated.toList()

  output:
  file "*blood_tissue_filtered.tsv" into blood_tissue_filtered mode flatten

  shell:
  if(params.blood_tissue_filter){
    '''
    blood_tissue_filter.R
    '''
  }
  '''
  for file in *.tsv
  do
    ln -s "$file" "${file/.tsv/_blood_tissue_filtered.tsv}"
  done
  '''
}

process annotation {

  tag { SM_tag }

  input:
  file filt from blood_tissue_filtered

  output:
  file "*hg38_multianno.txt" into annotated

  shell:
  SM_tag = filt.baseName.substring(0,12)
  '''
  table_annovar.pl -nastring NA -buildver hg38 --thread !{params.annovar_cpu} --onetranscript -remove -protocol refGene,exac03nontcga,esp6500siv2_all,1000g2015aug_all,gnomad_exome,clinvar_20170905,revel,cosmic80,avsnp147,ljb26_all -operation g,f,f,f,f,f,f,f,f,f -otherinfo !{filt} !{params.annovar_db}
  sed -i '1s/Otherinfo/QUAL\tFILTER\tINFO\tFORMAT\tGT\tIndividual\tStudy/' !{filt}.hg38_multianno.txt
  '''

}

process merge {

  publishDir params.out_folder, mode: 'move'

  input:
  file annotated_table from annotated.toList()

  output:
  file "*.tsv" into reformated_for_annovar mode flatten

  shell:
  '''
  mlr --tsv cat *.txt > big.tsv
  study=`head -n2 big.tsv | awk '{print $NF}' | grep -v "Study"`
  mv big.tsv TCGA_platypus_reformat_${study}.tsv
  '''
}

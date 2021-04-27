def helpMessage() {
    log.info"""
    ============================================================================
     :  Git version: ${version}
    ============================================================================
    Usage:

       Mandatory arguments:
           --out_dir        Output Directory (for things that aren't intermediate files to keep)
           --vcf            VCF to annotate w/snpEff
    """.stripIndent()
}

//paramater parsing and error checking:

workingDir = new File(".").getCanonicalPath()

baseDir = file("$baseDir")

def proc_git = "git -C $baseDir rev-parse HEAD".execute()
version = proc_git.text.trim()

params.help = false
if (params.help){
    println(params)
    helpMessage()
    exit 0
}

vcf_input = Channel.fromPath(params."vcf")
              .map {it -> [it.baseName.replaceFirst(/.vcf.gz$/, ""), it]}

process snpeff_annotate {
    publishDir params."out_dir" + '/snps-snpEff/all/', mode: 'copy', overwrite: true, pattern: '*.snpEff.ann.vcf.gz*'
    publishDir params."out_dir" + '/snps-snpEff/HIGH/', mode: 'copy', overwrite: true, pattern: '*.snpEff.ann.HIGH-IMPACT.vcf.gz*'
    publishDir params."out_dir" + '/snps-snpEff/MODERATE/', mode: 'copy', overwrite: true, pattern: '*.snpEff.ann.MODERATE-IMPACT.vcf.gz*'
    publishDir params."out_dir" + '/snps-snpEff/LOW/', mode: 'copy', overwrite: true, pattern: '*.snpEff.ann.LOW-IMPACT.vcf.gz*'
    publishDir params."out_dir" + '/snps-snpEff/MODIFIER/', mode: 'copy', overwrite: true, pattern: '*.snpEff.ann.MODIFIER-IMPACT.vcf.gz*'

    publishDir params."out_dir" + '/snps-snpEff/all/', mode: 'copy', overwrite: true, pattern: '*csv*'
    publishDir params."out_dir" + '/snps-snpEff/all/', mode: 'copy', overwrite: true, pattern: '*html*'


    container 'medicinalgenomics/snpeff-psicub-envy'
    memory '8 G'
    cpus 1

    input:
 
    set val(label), file(vcf_input) from vcf_input

    output:
    set file("*.snpEff.ann.vcf.gz"), file("*.snpEff.ann.vcf.gz.tbi"), file("*.csv"), file("*.summary.html"), file("*HIGH*"), file("*MODERATE*"), file("*LOW*"), file("*MODIFIER*") into snpeff_annotate_out

    script:
    memory = task.memory.toGiga()
    """
    sleep 0
    echo "Psicub_Envy.genome: Psicub_Envy" >> /opt/snpEff/snpEff.config

    java -Xmx${memory}g \
     -jar \$SNPEFF_HOME/snpEff.jar \
     Psicub_Envy \
     -csvStats ${label}.snpEff.csv \
     -nodownload \
     -canon \
     -v $vcf_input | bgzip -c > ${label}.snpEff.ann.vcf.gz

    tabix ${label}.snpEff.ann.vcf.gz

    (tabix -H ${label}.snpEff.ann.vcf.gz; gunzip -c ${label}.snpEff.ann.vcf.gz | grep \"|HIGH|\") | bgzip -c > ${label}.snpEff.ann.HIGH-IMPACT.vcf.gz
    tabix ${label}.snpEff.ann.HIGH-IMPACT.vcf.gz

    (tabix -H ${label}.snpEff.ann.vcf.gz; gunzip -c ${label}.snpEff.ann.vcf.gz | grep \"|MODERATE|\") | bgzip -c > ${label}.snpEff.ann.MODERATE-IMPACT.vcf.gz
    tabix ${label}.snpEff.ann.MODERATE-IMPACT.vcf.gz

    (tabix -H ${label}.snpEff.ann.vcf.gz; gunzip -c ${label}.snpEff.ann.vcf.gz | grep \"|LOW|\") | bgzip -c > ${label}.snpEff.ann.LOW-IMPACT.vcf.gz
    tabix ${label}.snpEff.ann.LOW-IMPACT.vcf.gz

    (tabix -H ${label}.snpEff.ann.vcf.gz; gunzip -c ${label}.snpEff.ann.vcf.gz | grep \"|MODIFIER|\") | bgzip -c > ${label}.snpEff.ann.MODIFIER-IMPACT.vcf.gz
    tabix ${label}.snpEff.ann.MODIFIER-IMPACT.vcf.gz    

    mv snpEff_summary.html ${label}.snpEff.summary.html
    """

}

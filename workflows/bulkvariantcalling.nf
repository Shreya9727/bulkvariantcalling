/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { BWA_INDEX } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM } from '../modules/nf-core/bwa/mem/main'
include { GATK4_MARKDUPLICATES } from '../modules/nf-core/gatk4/markduplicates/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary'
include { GATK4_HAPLOTYPECALLER } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_bulkvariantcalling_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BULKVARIANTCALLING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main: 
    def ch_versions      = channel.empty()
    def ch_multiqc_files = channel.empty()

    ch_fasta       = Channel.value([[id: 'reference'], file(params.fasta)])
    ch_fasta_plain = Channel.value(file(params.fasta))
    ch_fasta_meta  = Channel.value([[id: 'genome'], file(params.fasta)])

    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    BWA_INDEX(ch_fasta)
    BWA_MEM(
        ch_samplesheet,
        BWA_INDEX.out.index,
        ch_fasta,
        true
	)

    SAMTOOLS_FAIDX(
	    ch_fasta_meta.map { meta, fasta -> [ meta, fasta, [] ] },
	    false
    )

    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta_meta)   

   
    ch_fai_meta = SAMTOOLS_FAIDX.out.fai.map { meta, fai -> [[id: 'genome'], fai] }.first()
    ch_dict_meta = GATK4_CREATESEQUENCEDICTIONARY.out.dict.map { meta, dict -> [[id: 'genome'], dict] }.first()


    GATK4_MARKDUPLICATES(
       BWA_MEM.out.bam,
       ch_fasta_plain,
       SAMTOOLS_FAIDX.out.fai.map { meta, fai -> fai }.first()
	)

    ch_markdup_bam_bai = GATK4_MARKDUPLICATES.out.bam
	    .join(GATK4_MARKDUPLICATES.out.bai)
	    .map { meta, bam, bai -> [ meta, bam, bai, [], [] ] }   
    GATK4_HAPLOTYPECALLER(
	    ch_markdup_bam_bai,
	    ch_fasta_meta,
	    ch_fai_meta,
    	    ch_dict_meta,
    	    [[:], []],   // dbsnp — empty/optional for now
            [[:], []]    // dbsnp_tbi — empty/optional for now
	)
    
    MULTIQC(
	    FASTQC.out.zip.mix(GATK4_MARKDUPLICATES.out.metrics)
	        .map { meta, files -> files }   // drop meta, keep just the file(s)
	        .flatten()                       
	        .collect()
	        .map { files ->
	            [
	                [id: 'bulkvariantcalling'],
	                files,
	                multiqc_config
	                    ? file(multiqc_config, checkIfExists: true)
	                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
	                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
	                [],
	                [],
	            ]
	        }
	)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'bulkvariantcalling_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

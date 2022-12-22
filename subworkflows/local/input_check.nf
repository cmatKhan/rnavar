//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    ch_reads = Channel.empty()

    SAMPLESHEET_CHECK (
        samplesheet
    )
    .csv
    .splitCsv ( header:true, sep:',' )
    .multiMap{ it ->
        def parsed_row = create_input_channel(it)
        reads: parsed_row.fastq
        bam_bai: parsed_row.bam_bai
    }
    .set { ch_parsed_input }

    emit:
    reads = ch_parsed_input.reads             // channel: [ val(meta), [ reads ] ]
    bam_bai = ch_parsed_input.bam_bai         // channel: [ val(meta), file(bam), file(bai) ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_input_channel(LinkedHashMap row) {
    // create output map
    def output_map = [:]
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    // note: do not set single_end flag here -- instead, set it in the fastq case
    // below
    // create fastq and bam_bai-- this is for the case in which there are
    // bams, but no fastqs
    output_map.bam_bai = [meta, null, null]
    output_map.fastq = [meta,[]]

    // check some assumptions on the samplesheet, namely that each row must have
    // at minimum a valid fastq_1 OR valid bam, but not both
    if(row.fastq_1 == ''){
        if(!file(row.bam).exists()){
            exit 1, "ERROR: Please check input samplesheet -> fastq_1 is blank and bam file does not exist:\n${row.bam}"
        } else {
            if (!file(row.bam+'.bai').exists()){
                exit 1, "ERROR: Please check input samplesheet -> bam must have a .bai index file in same dir. One does not exist:\n${row.bam+'.bai'}"
            } else {
                output_map.bam_bai = [ meta, file(row.bam), file(row.bam+'.bai') ]
            }
        }
    } else if(row.bam == ''){

        // set single end boolean in fastq case
        meta.single_end = row.single_end.toBoolean()

        if(!file(row.fastq_1).exists()){
            exit 1, "ERROR: Please check input samplesheet -> bam is blank and Read 1 does not exist:\n ${row.bam}"
        } else if (meta.single_end){
                output_map.fastq = [ meta, [ file(row.fastq_1) ] ]
        } else {
            if(file(row.fastq_2).exists()){
                output_map.fastq = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
            } else {
                exit 1, "ERROR: Please check input samplesheet -> single_end is set to false but Read 2 does not exist:\n${row.fastq_1}"
            }
        }
    }

    return output_map
}

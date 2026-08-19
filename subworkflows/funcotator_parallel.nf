include { split_funcotator_chunks } from '../modules/local/annotation/small_variants/funcotator/split_funcotator_chunks/main.nf'
include { run_funcotator_parallel } from '../modules/local/annotation/small_variants/funcotator/run_funcotator_parallel/main.nf'
include { run_funcotator_single } from '../modules/local/annotation/small_variants/funcotator/run_funcotator_single/main.nf'
include { merge_funcotator_chunks } from '../modules/local/annotation/small_variants/funcotator/merge_funcotator_chunks/main.nf'

workflow FUNCOTATOR {

    take:
        ch_input
        safe_funcotator_germline_db
        safe_funcotator_somatic_db
        safe_transcript_list
        funcotator_n_chunks

    main:

        if( funcotator_n_chunks > 1 ) {

            split_flat_ch = split_funcotator_chunks(ch_input, funcotator_n_chunks)
                .flatMap { meta, chunk_index, n_parts_file, vcfs, vcf_full, maf_skip ->

                    def key = groupKey(
                        "${meta.patient}_${chunk_index}",
                        n_parts_file.text.trim().toInteger()
                    )

                    (vcfs instanceof List ? vcfs : [vcfs]).collect { vcf_file ->
                        [
                            key,
                            meta,
                            chunk_index,
                            (vcf_file.simpleName =~ /part_(\d+)/)[0][1],
                            vcf_file,
                            vcf_full,
                            maf_skip
                        ]
                    }
                }

            grouped_ch = run_funcotator_parallel(
                    split_flat_ch,
                    safe_funcotator_germline_db,
                    safe_funcotator_somatic_db,
                    safe_transcript_list
                )
                .groupTuple(by: 0)
                .map { group_key, metas, chunk_indexes, mini_indexes, mafs, vcfs, vcf_fulls, maf_skips ->

                    [
                        metas[0],
                        chunk_indexes[0],
                        mafs,
                        vcfs,
                        vcf_fulls[0],
                        maf_skips[0]
                    ]
                }

            merged_ch = merge_funcotator_chunks(grouped_ch)

        } else {

            merged_ch = run_funcotator_single(
                ch_input,
                safe_funcotator_germline_db,
                safe_funcotator_somatic_db,
                safe_transcript_list
            )

        }

    emit:
        merged_ch
}

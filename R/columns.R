## Columns


# Objects
# # x_read_output_cols
# # escape_gene_cols


x_read_output_cols = c(
    'mouse_id',
    # 'chromosome_mat',  # always X in this format
    # 'chromosome_pat',  # always X in this format
    # 'locus_mat',  # came with Berletch's data
    # 'locus_pat',  # came with Berletch's data
    'num_total_reads',
    'bias_xi_div_xa',
    'gene_id_mat',
    'gene_id_pat',
    'exon_length_mat',  # required now
    'exon_length_pat',  # required now
    'gene_name',
    'num_reads_mat',
    'num_reads_pat',
    'ratio_xi_over_xa',
    # 'total_reads',  # this was an intermediate for calculation only
    'srpm_mat',
    'srpm_pat',
    # 'rpm_mat',  # don't need this
    # 'rpm_pat',  # don't need this
    'rpkm_mat',  # not present if RPKMs are merged
    'rpkm_pat',  # not present if RPKMs are merged
    'rpkm',
    'pct_xi',
    'corrected_pct_xi',
    'lower_confidence_interval',
    'upper_confidence_interval',
    'xi_srpm_gte_2',
    'rpkm_gt_1',
    'lower_ci_gt_0'
)


# shortlist for easy comparison
escape_gene_cols = c(
    'mouse_id',
    'gene_name',
    'num_reads_mat',
    'num_reads_pat',
    'ratio_xi_over_xa',
    'lower_confidence_interval',
    'upper_confidence_interval',
    'srpm_mat',
    'srpm_pat',
    'rpkm'
)

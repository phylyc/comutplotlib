class MutationAnnotation(object):
    """Defines all mutation annotation column names and entry names."""

    # columns to be added to the mutation annotation data
    channel = "Mutation_Channel"
    cropped_context = "Cropped_Context"

    signature = "Signature"

    # mutation types
    snv = "SNP"  # single nucleotide polymorphism (variant)
    dnv = "DNP"  # double nucleotide polymorphism (variant)
    tnv = "TNP"  # triple nucleotide polymorphism (variant)
    mnv = "MNP"  # multi nucleotide polymorphism (variant)
    insertion = "INS"
    deletion = "DEL"

    # mutation effects
    # https://gatk.broadinstitute.org/hc/en-us/articles/360037224432-Funcotator
    # Variant lies between exons within the bounds of the chosen transcript. Only valid
    # for Introns.
    intron = "Intron"
    # Variant is on the 5'UTR for the chosen transcript. Only valid for UTRs.
    utr5 = "5'UTR"
    # Variant is on the 3'UTR for the chosen transcript. Only valid for UTRs.
    utr3 = "3'UTR"
    # Intergenic region. Does not overlap any transcript. Only valid for IGRs.
    igr = "IGR"
    # The variant is upstream of the chosen transcript. Only valid for IGRs.
    # Variants within (by default) 5000 bases of the 5' end of a transcript (and not
    # overlapping any part of the transcript itself) will be annotated as being in the
    # 5' flanking region of that transcript.
    flank5 = "5'Flank"
    # The variant is downstream of the chosen transcript. Only valid for IGRs.
    # by default: 0 bases.
    flank3 = "3'Flank"
    # The point mutation alters the protein structure by one amino acid. Can occur in
    # Coding regions or Introns.
    missense = "Missense_Mutation"
    # A premature stop codon is created by the variant. Can occur in Coding regions or
    # Introns.
    nonsense = "Nonsense_Mutation"
    # Variant removes stop codon. Can occur in Coding regions or Introns.
    nonstop = "Nonstop_Mutation"
    # Variant is in coding region of the chosen transcript, but protein structure is
    # identical. Can occur in Coding regions or Introns.
    silent = "Silent"
    # The variant is within a configurable number of bases of a splice site. See the
    # secondary classification to determine if it lies on the exon or intron side. Can
    # occur in Coding regions or Introns.
    splice_site = "Splice_Site"
    # Deletion that keeps the sequence in frame. Can occur in Coding regions or Introns.
    in_frame_del = "In_Frame_Del"
    # Insertion that keeps the sequence in frame. Can occur in Coding regions or
    # Introns.
    in_frame_ins = "In_Frame_Ins"
    # Deletion that moves the sequence out of frame. Can occur in Coding regions or
    # Introns.
    frame_shift_del = "Frame_Shift_Del"
    # Insertion that moves the coding sequence out of frame. Can occur in Coding regions
    # or Introns.
    frame_shift_ins = "Frame_Shift_Ins"
    # Point mutation that overlaps the start codon. Can occur in Coding regions or
    # Introns.
    start_codon_snp = "Start_Codon_SNP"
    # Deletion that overlaps the start codon. Can occur in Coding regions or Introns.
    start_codon_del = "Start_Codon_Del"
    # Insertion that overlaps the start codon. Can occur in Coding regions or Introns.
    start_codon_ins = "Start_Codon_Ins"
    # New start codon is created by the given variant using the chosen transcript.
    # However, it is in frame relative to the coded protein, meaning that if the coding
    # sequence were extended then the new start codon would be in frame with the
    # existing start and stop codons. This can only occur in a 5' UTR.
    de_novo_start_in_frame = "De_novo_Start_InFrame"
    # New start codon is created by the given variant using the chosen transcript.
    # However, it is out of frame relative to the coded protein, meaning that if the
    # coding sequence were extended then the new start codon would NOT be in frame with
    # the existing start and stop codons. This can only occur in a 5' UTR.
    de_novo_start_out_of_frame = "De_novo_Start_OutOfFrame"
    # Variant lies on one of the RNA transcripts. (special catch-all case)
    rna = "RNA"
    # Variant lies on one of the lincRNAs. (special catch-all case)
    linc_rna = "lincRNA"
    read_through = "Read-through"
    stop_codon_ins = "Stop_Codon_Ins"
    stop_codon_del = "Stop_Codon_Del"
    translation_start_site = "Translation_Start_Site"

    # pooled effects
    structural = "Structural"
    synonymous = "Synonymous"
    nonsynonymous = "Non-synonymous"
    gain_of_function = "Gain of Function"
    loss_of_function = "Loss of Function"

    # strand
    plus_strand = "+"
    minus_strand = "-"
    both_strands = "+/-"

    # columns in the mutation annotation file
    patient = "Patient"
    sample = "SIF_sample_name"
    sample_barcode = "Tumor_Sample_Barcode"
    gene_name = "Hugo_Symbol"
    gene_id = "HGNC_Ensembl_Gene_ID"
    # entrez_id = "Entrez_Gene_Id"
    annotation_transcript = "Annotation_Transcript"
    chromosome = "Chromosome"
    start_pos = "Start_Position"
    end_pos = "End_Position"
    strand = "Transcript_Strand"
    transcript_exon = "Transcript_Exon"
    transcript_pos = "Transcript_Position"
    type = "Variant_Type"
    effect = "Variant_Classification"
    mutation_status = "Mutation_Status"
    ref_allele = "Reference_Allele"  # reference allele
    alt_allele = "Tumor_Seq_Allele2"  # alternate allele
    matched_normal_ref_allele = "Match_Norm_Seq_Allele1"
    matched_normal_alt_allele = "Match_Norm_Seq_Allele2"
    genome_change = "Genome_Change"
    codon_change = "Codon_Change"
    protein_change = "Protein_Change"
    tumor_fraction = "tumor_f"
    ref_count = "t_ref_count"  # number of reads of the reference allele
    alt_count = "t_alt_count"  # number of reads of the alternate allele
    matched_normal_ref_count = "n_ref_count"  # number of reads of the reference allele
    matched_normal_alt_count = "n_alt_count"  # number of reads of the alternate allele
    conservation = "Conservation"
    context = "ref_context"  # by default the window size is 10 bases.
    # The window size does not include any bases in the variant alternate allele itself.
    # By default, the window size is 200 bases.
    gc_content = "gc_content"

    clin_var_filter = "ClinVar_VCF_CLNSIG"
    dbSNP_filter = "dbSNP_FILTER"
    gnomad_exome_filter = "gnomAD_exome_FILTER"
    gnomad_genome_filter = "gnomAD_genome_FILTER"
    gnomad_exome_AF = "gnomAD_exome_AF"
    gnomad_genome_AF = "gnomAD_genome_AF"
    gnomad_exome_AF_popmax = "gnomAD_exome_AF_popmax"
    gnomad_genome_AF_popmax = "gnomAD_genome_AF_popmax"

    # if splice_site is intron or exon
    # Gencode version depends on annotation source!
    splice_site_class = "Gencode_34_secondaryVariantClassification"
    # Difference in alignment score between best and next-best alignment
    best_alignment_difference = "ALIGN_DIFF"
    # "AD" is the count of informative reads that support a given haplotype.
    # This does not count reads at the site that were uninformative to any of
    # the alleles. This means that reads that are ambiguous, or don't span the
    # entire repeat at STR sites will likely be excluded from AD, but might
    # still be in the DP. This can often lead to a disagreement with DP.
    allelic_depth = "AD"
    # approximate read depth
    # refers to the number of reads that were set to the genotyper at a given site
    read_depth = "DP"
    # ???
    # "ECNT"
    # Phred-scaled quality that alt alleles are not germline variants
    not_germline_quality = "GERMQ"
    # median base quality by allele
    median_base_quality_by_allele = "MBQ"
    # Median fragment length of reads supporting each allele
    median_fragment_length_by_allele = "MFRL"
    # median mapping quality by allele
    median_mappability_quality_by_allele = "MMQ"
    # median distance from end of read
    median_end_of_read_distance = "MPOS"
    # Number of joint alignments
    number_of_joint_alignments = "NALIGNS"
    # Negative log 10 odds of artifact in normal with same allele fraction as tumor
    nlog10_odds_is_artifact_in_normal = "NALOD"
    # Normal log 10 likelihood ratio of diploid het versus hom alt genotype
    nlog10_diploid = "NLOD"
    # negative log 10 population allele frequency of alt alleles
    nlog10_population_alt_allele_frequency = "POPAF"
    # Phred-scaled qualities that alt allele are not due to read orientation artifact
    read_orientation_quality = "ROQ"
    # number of times tandem repeat unit is repeated, for each allele (including ref)
    repeat_size_by_allele = "RPA"
    # tandem repeat unit (bases)
    repeat_unit = "RU"
    # variant is a short tandem repeat
    is_tandem_repeat = "STR"
    # Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors
    slippage_quality = "STRQ"
    # log 10 likelihood ratio of variant existing versus not existing
    log10_odds_exists = "TLOD"
    # ???
    # "UNITIGS"

    default_columns = [
        patient,
        sample,
        gene_name,
        gene_id,
        chromosome,
        start_pos,
        end_pos,
        type,
        effect,
        ref_allele,
        alt_allele,
        context,
    ]

    column_dtype = {
        # patient: str,
        # sample: str,
        # gene_name: str,
        # gene_id: str,
        chromosome: str,  # category
        # start_pos: int,
        # end_pos: int,
        # effect: str,  # category
        # ref_allele: str,  # category
        # alt_allele: str,  # category
        # context: str,
        # ref_count: int,
        # alt_count: int,
    }

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

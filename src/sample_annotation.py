class SampleAnnotation(object):

    sample = "sample_id"
    patient = "clinical_id"
    platform = "platform"
    platform_abv = "Platform"
    data_type = "data_type"
    center = "center"
    bam_path = "clean_bam_file_capture"
    gc_path = "google_cloud_path"
    cancer_type = "primary_cancer_type"
    histotype = "primary_histotype"
    histology = "Histology"
    sample_type_long = "sample_type"
    sample_type = "sample_type_abv"
    sample_description = "sample_description"
    material = "specimen_material"
    protocol = "Protocol(s)"
    sex = "sex"
    sex_genotype = "SEX_GENOTYPE"

    contamination = "contamination"
    contamination_error = "contamination_error"
    tumor_purity = "tumor_purity"
    tumor_purity_error = "tumor_purity_error"
    tmb = "tmb"
    tmb_error = "tmb_error"

    age_p_dx = "age_p_dx"
    age_bm_dx = "age_bm_dx"
    age_death = "age_death"
    months_p_dx_to_bm_dx = "months_p_dx_to_bm_dx"
    months_p_dx_to_death = "months_p_dx_to_death"
    months_bm_dx_to_death = "months_bm_dx_to_death"
    date_p_rx = "date_p_rx"
    date_bm_rx = "date_bm_rx"
    smoking_hx = "smoking_hx"

    whole_exome_seq = "WES"
    whole_genome_seq = "WGS"

    ultra_low_pass = "ULP"
    rna = "RNA"

    default_columns = [patient, sample]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

from adjustText import adjust_text
from collections import Counter, defaultdict
from itertools import chain
import matplotlib.pyplot as plt
from matplotlib import colors, patches, ticker
import numpy as np
import pandas as pd
import re

from src.amino_acid import AminoAcidChange
from src.functional_effect import sort_functional_effects
from src.math import decompose_rectangle_into_polygons
from src.mutation_annotation import MutationAnnotation as MutA
from src.plotter import Plotter
from src.sample_annotation import SampleAnnotation


def plot_comutation_table(
    *args,
    output: str = ".",
    maf=None,
    pool_as=None,
    variant_filter=None,
    sif=None,
    meta_data_rows=None,
    meta_data_rows_per_sample=("Sample Type", "Material", "Platform", "Contamination", "Tumor Purity"),
    gistic=None,
    by: str = MutA.patient,
    label_columns: bool = True,
    column_order: tuple = None,
    index_order: tuple = None,
    snv_comut_threshold: int = None,
    cnv_comut_threshold: int = None,
    snv_quantile: float = 1 - 50 / 20000,
    cnv_quantile: float = 1 - 10 / 20000,
    select_quantile: bool = False,
    interesting_genes: set = None,
    snv_interesting_genes: set = None,
    cnv_interesting_genes: set = None,
    interesting_gene: str = None,
    interesting_gene_comut_threshold_percent: float = None,
    ground_truth_genes: dict[str, list[str]] = None,
    low_amp_threshold: int = 1,
    high_amp_threshold: int = 2,
    low_del_threshold: int = -1,
    high_del_threshold: int = -2,
    total_prevalence_threshold: float = None,
    model_names: list = ("SNV burden model", "CNV burden model"),
    palette: dict[str, tuple[float]] = None,
    parts_to_plot: list = None,
    max_xfigsize: int = None,
    yfigsize: int = None,
    sub_folders=(),
    file_name_prefix="",
    file_name_suffix="",
    **kwargs
):
    return None

    plotter = Plotter(out_dir=output)

    _palette = plotter.palette.copy()
    if palette is not None:
        # overwrite native palette with specified palette:
        for key, color in palette.items():
            plotter.palette[key] = color
    palette = _palette

    if meta_data_rows is None:
        meta_data_rows = [
            "Sample Type",
            "Material",
            "Contamination",
            "Tumor Purity",
            "Platform",
            "has matched N",
            "Sex",
            "Histology",
        ]
    if parts_to_plot is None:
        parts_to_plot = [
            "mutation burden",
            "mutation burden legend",
            "recurrence",
            "prevalence",
            "total prevalence",
            "cytoband",
            "comutations",
            "comutations legend",
            "meta data",
            "meta data legend",
        ]
    ######################
    # PROCESSING OF DATA #
    ######################

    if by == MutA.sample:
        columns = sif.samples
        columns_name = SampleAnnotation.sample
    elif by == MutA.patient:
        columns = sif.patients
        columns_name = SampleAnnotation.patient
    else:
        raise ValueError("The argument 'by' needs to be one of {" + f"{MutA.sample}, {MutA.patient}" + "} but received " + f"{by}.")

    uninteresting_effects = [
        # MutA.utr5, MutA.utr3, MutA.flank5, MutA.flank3,
        MutA.synonymous, MutA.silent, MutA.igr, MutA.intron
    ]

    maf = (
        maf
        .select(selection={by: list(columns)})
        .select(
            selection={MutA.effect: uninteresting_effects},
            complement=True
        ).select_variants(
            variant_filter=variant_filter,
            skip=[
                maf.variant_filter.num_variants_per_patient_filter,
                maf.variant_filter.vaf_filter,
            ],
            verbose=False
        )
    )

    base_protein_change_col = "__base_protein_change"
    maf.assign_column(
        name=base_protein_change_col,
        value=maf.data[MutA.protein_change].apply(
            lambda p: (lambda aac: aac.base.code + aac.pos)(AminoAcidChange.from_string(p))
            if isinstance(p, str)
            else pd.NA
        ),
        inplace=True,
    )

    maf_patient_dedup = maf.drop_patient_duplicates()

    maf_pos_aggregate = maf_patient_dedup.stratify_mutation_counts(
        by=[MutA.gene_name, base_protein_change_col],
        expand_index=False,
        counts_column_name=None,
    )

    pre_cna = gistic.sample_table.reindex(columns=columns)
    has_high_cnv = ~pre_cna.isna() & (pre_cna.ge(high_amp_threshold) | pre_cna.le(high_del_threshold))

    ########################
    # SELECT GENES TO PLOT #
    ########################

    def get_interesting_genes(aggregate, threshold=None, func=min, q=0.999):
        threshold = (
            int(
                func(
                    np.quantile(aggregate.groupby(MutA.gene_name).agg("max"), q=q),
                    max(0.03 * len(columns), 2),
                )
            )
            if threshold is None
            else threshold
        )
        _interesting_genes = (
            aggregate.loc[aggregate >= threshold]
            .index.get_level_values(level=MutA.gene_name)
            .unique()
        )
        return threshold, set(_interesting_genes)

    if interesting_genes is None:
        if snv_interesting_genes is None and interesting_gene is None:
            snv_comut_threshold, snv_interesting_genes = get_interesting_genes(
                maf_pos_aggregate, snv_comut_threshold, min, snv_quantile
            )
        if cnv_interesting_genes is None and interesting_gene is None:
            cnv_comut_threshold, cnv_interesting_genes = get_interesting_genes(
                has_high_cnv.sum(axis=1), cnv_comut_threshold, max, cnv_quantile
            )
        if snv_interesting_genes is None and cnv_interesting_genes is None and interesting_gene is not None:
            has_snv = (
                maf
                .data
                .groupby(by=[by, MutA.gene_name])
                .agg("size")
                .astype(bool)
                .unstack(by)
                .reindex(index=pre_cna.index)
                .reindex(columns=columns)
                .fillna(False)
            )
            has_mut = has_snv | has_high_cnv
            gene = has_mut.loc[interesting_gene]
            good_genes = has_mut.apply(
                lambda g: (g & gene).sum() / gene.sum() >= interesting_gene_comut_threshold_percent,
                axis=1
            )
            snv_interesting_genes = set(good_genes.loc[good_genes].index)
            cnv_interesting_genes = snv_interesting_genes
            snv_comut_threshold = int(gene.sum() * interesting_gene_comut_threshold_percent) + 1
            cnv_comut_threshold = int(gene.sum() * interesting_gene_comut_threshold_percent) + 1

        interesting_genes = snv_interesting_genes | cnv_interesting_genes

    if ground_truth_genes is not None:
        for _, ground_truth_gene_list in ground_truth_genes.items():
            interesting_genes |= set(ground_truth_gene_list)

    if len(interesting_genes) == 0:
        return None

    ###############################
    # CREATING DATAFRAMES TO PLOT #
    ###############################

    mut_burden_data = (
        (maf_patient_dedup if by == MutA.patient else maf)
        .pool_annotations(pool_as=pool_as)
        .stratify_mutation_counts(by=[by, MutA.effect], counts_column_name=None)
        .sparse.to_dense()
        .unstack(level=MutA.effect)
        .reindex(index=columns)
        .fillna(0)
    )
    mut_recurrence = (
        maf_patient_dedup
        .data[[MutA.gene_name, MutA.effect, MutA.protein_change]]
        .replace(np.nan, "NaN")
        .groupby([MutA.gene_name, MutA.effect, MutA.protein_change])
        .agg("size")
    )
    mut_prevalence_counter = (
        maf_patient_dedup
        .select({MutA.gene_name: list(interesting_genes)})
        .data
        .drop_duplicates(subset=[MutA.patient, MutA.gene_name, MutA.effect])
        .groupby(by=[MutA.gene_name])
        .agg({MutA.effect: Counter})
        .get(MutA.effect)
        .reindex(index=interesting_genes)
    )
    mut = (
        (maf_patient_dedup if by == MutA.patient else maf)
        .select({MutA.gene_name: list(interesting_genes)})
        .data
        .groupby(by=[by, MutA.gene_name])
        .agg({MutA.effect: sort_functional_effects})
        .get(MutA.effect)
        .unstack(by)
        .reindex(index=interesting_genes)
        .reindex(columns=columns)
    )
    cna = (
        gistic
        .sample_table
        .reindex(index=interesting_genes)
        .reindex(columns=columns)
    )

    ##############################
    # FURTHER SELECTION OF GENES #
    ##############################

    if select_quantile:
        has_snv = ~mut.isna()
        has_high_cnv = ~cna.isna() & (cna.ge(high_amp_threshold) | cna.le(high_del_threshold))

        gene_counts = (has_snv.astype(int) + has_high_cnv.astype(int)).sum(axis=1)
        gene_count_threshold = np.quantile(gene_counts, 1 / 3)
        genes = gene_counts.loc[gene_counts > gene_count_threshold].index
    else:
        genes = list(interesting_genes)

    mut = mut.loc[genes]
    cna = cna.loc[genes]
    cytobands = gistic.data["Cytoband"].reindex(genes, axis=0).fillna("X")

    if mut.empty & cna.empty:
        return None

    ############################
    # SORTING ROWS AND COLUMNS #
    ############################

    has_snv = ~mut.isna()
    has_high_amp = ~cna.isna() & cna.ge(high_amp_threshold)
    has_high_del = ~cna.isna() & cna.le(high_del_threshold)
    has_low_amp = ~cna.isna() & cna.eq(low_amp_threshold)
    has_low_del = ~cna.isna() & cna.eq(low_del_threshold)
    has_high_cnv = has_high_amp | has_high_del
    has_low_cnv = has_low_amp | has_low_del
    has_mut = has_snv | has_high_cnv

    # sort genes
    sorted_features = (
        (has_snv.astype(int) + has_high_cnv.astype(int))
        .sum(axis=1)
        .to_frame("mut_count")
        .join(has_snv.any(axis=1).to_frame("has_snv"))
        .join(has_low_cnv.astype(int).sum(axis=1).to_frame("low_cnv_count"))
        .join(has_low_cnv.any(axis=1).to_frame("has_low_cnv"))
        .sort_values(by=["mut_count", "has_snv", "has_low_cnv", "low_cnv_count"], ascending=True)
    )
    genes = sorted_features.index.to_list()
    # bring the interesting gene to the front if the heatmap:
    if interesting_gene is not None:
        genes.pop(genes.index(interesting_gene))
        genes += [interesting_gene]

    # group genes in the same cytoband together since they are co-amplified or co-deleted.
    cytobands = cytobands.loc[genes]
    cytoband_groups = cytobands[::-1].drop_duplicates().values
    cytoband_key = {cb: i for i, cb in enumerate(cytoband_groups)}

    gene_key = {}
    feature_count = 0
    previous_row_tuple = None
    for gene, row in sorted_features.loc[genes[::-1]].iterrows():
        row_tuple = tuple(row)
        if previous_row_tuple is None or row_tuple != previous_row_tuple:
            feature_count += 1
        gene_key[gene] = feature_count
        previous_row_tuple = row_tuple

    genes = cytobands.reset_index().set_axis(genes, axis=0).apply(
        lambda row: (cytoband_key[row["Cytoband"]], gene_key[row[MutA.gene_name]], row[MutA.gene_name]),
        axis=1
    ).sort_values(ascending=False).index.to_list()

    # bring the interesting gene to the front (again):
    # and remove genes in the same or neighboring cytobands
    if interesting_gene is not None:
        # genes.pop(genes.index(interesting_gene))
        sorted_cytoband_groups = sorted(
            cytoband_groups,
            key=lambda cytoband: [c if i % 2 else int(c) for i, c in enumerate(re.split('(\d+)', cytoband)[1:-1])]
        )
        sorted_cytoband_key = {cb: i for i, cb in enumerate(sorted_cytoband_groups)}
        cytoband_diff = cytobands.apply(lambda c: sorted_cytoband_key[c]) - sorted_cytoband_key[cytobands.loc[interesting_gene]]
        close_genes = cytoband_diff.abs() < 5
        for g in close_genes.loc[close_genes].index:
            genes.pop(genes.index(g))
        genes += [interesting_gene]

    # sort columns
    has_burden = mut_burden_data.gt(0).any(axis=1).astype(int).to_frame("burden").T
    has_cnv = (cna.ne(0) & ~cna.isna()).any(axis=0).astype(int).to_frame("has_cnv").T
    has_high_mut = (
        4 * has_snv.loc[genes].astype(int)
        + 3 * (has_snv & ~has_high_cnv).loc[genes].astype(int)
        + 2 * has_high_amp.loc[genes].astype(int)
        + has_high_del.loc[genes].astype(int)
    )
    has_low_mut = (
        + 2 * has_low_amp.loc[genes].astype(int)
        + has_low_del.loc[genes].astype(int)
    )
    columns = (
        pd.concat([has_low_mut, has_cnv, has_burden, has_high_mut]).T
        .apply(lambda x: tuple(reversed(tuple(x))), axis=1)
        .sort_values(ascending=False)
        .index
    ) if column_order is None else [c for c in column_order if c in columns]

    total_prevalence_per_gene = has_mut.loc[genes].astype(int).sum(axis=1) / len(columns)

    # remove genes below a percentage of prevalence:
    if total_prevalence_threshold is not None:
        total_prevalence_per_gene = total_prevalence_per_gene.loc[
            total_prevalence_per_gene >= total_prevalence_threshold
        ]
        genes = total_prevalence_per_gene.index

    mut = mut.loc[genes, columns]
    cna = cna.loc[genes, columns]
    total_prevalence_overall = has_mut.any(axis=0).astype(int).sum() / len(columns)
    mut_prevalence_counter = mut_prevalence_counter.loc[genes]
    cna_counter = cna.apply(Counter, axis=1)
    cytobands = cytobands.loc[genes]
    if snv_interesting_genes is not None and cnv_interesting_genes is not None:
        model_annotation = pd.DataFrame(
            [[g in snv_interesting_genes, g in cnv_interesting_genes] for g in genes],
            index=genes,
            columns=["snv", "cnv"]
        )
    else:
        model_annotation = None
    mut_burden_data = mut_burden_data.loc[columns]
    # mut_protein_data = pd_util.subframe(
    #     df=mut_recurrence, batch=genes, batch_level=MutA.gene_name
    # )
    mut_protein_data = {
        g: (
            mut_recurrence.loc[g]
            if g in mut_recurrence.index.get_level_values(level=MutA.gene_name)
            else None
        )
        for g in genes
    }

    if mut.empty | cna.empty:
        return None

    ###################
    # SORTING EFFECTS #
    ###################

    effects = set(chain.from_iterable([l for _, l in np.ndenumerate(mut) if isinstance(l, list)]))
    effects = sort_functional_effects(effects)

    ##############
    # COLOR MAPS #
    ##############

    def better_effect_legend(effect):
        return (
            effect
            .replace("_", " ")
            .replace("OutOfFrame", "ooF")
            .replace("InFrame", "iF")
        )

    mut_burden_legend_cmap = {
        better_effect_legend(effect): palette.get(effect, palette.grey)
        for effect in reversed(mut_burden_data.columns)
    }
    # Some functional effects are painted with the same color, which messes up
    # the legend, so we have to construct the legend color map separately.
    mut_cmap = {}
    inv_mut_cmap = defaultdict(list)
    for effect in effects:
        color = palette.get(effect, palette.grey)
        mut_cmap[effect] = color
        inv_mut_cmap[color].append(better_effect_legend(effect))
    mut_legend_cmap = {
        " / ".join(e_list): c for c, e_list in inv_mut_cmap.items()
    }

    _cna_cmap = {
        high_amp_threshold: palette.red,
        low_amp_threshold: palette.adjust_lightness(palette.lightred, 1.16),
        0: palette.backgroundgrey,
        low_del_threshold: palette.adjust_lightness(palette.lightblue, 1.16),
        high_del_threshold: palette.blue,
    }
    _cna_names = [
        "High Amplification",
        "Low Amplification",
        "Baseline",
        "Low Deletion",
        "High Deletion",
    ]
    cna_cmap = {}
    cna_names = []
    for name, (key, color) in zip(_cna_names, _cna_cmap.items()):
        if key in np.unique(cna):
            cna_cmap[key] = color
            cna_names.append(name)

    #############
    # META DATA #
    #############

    # meta_columns = [
    #     SIF.sample,
    #     SIF.patient,
    #     SIF.platform_abv,
    #     SIF.sex,
    #     SIF.sample_type,
    #     SIF.material,
    #     SIF.histology,
    #     "Paired",
    #     "contamination",
    #     "tumor_purity",
    # ]
    # meta_columns = [c for c in meta_columns if c in sif.data.columns]
    # per_sample_columns = [
    #     c
    #     for c in [SIF.sample_type, SIF.material, SIF.platform_abv, "contamination", "tumor_purity"]
    #     if c in meta_columns
    # ]
    # per_patient_columns = [c for c in meta_columns if c not in per_sample_columns]
    # sif_meta_data = (
    #     sif.select({columns_name: list(columns)})
    #     .data[meta_columns]
    #     .groupby(by=[columns_name])
    #     .agg({k: list for k in per_sample_columns} | {k: lambda l: np.unique(l)[0] for k in per_patient_columns})
    #     .reindex(index=columns)
    # )
    sif_meta_data = sif.select({columns_name: list(columns)}).data

    sif_meta_data["Sample Type"] = sif_meta_data[SampleAnnotation.sample_type]
    sif_meta_data["Material"] = sif_meta_data[SampleAnnotation.material]
    sif_meta_data["Contamination"] = sif_meta_data["contamination"]
    sif_meta_data["Tumor Purity"] = sif_meta_data["tumor_purity"]
    sif_meta_data["Platform"] = sif_meta_data[SampleAnnotation.platform_abv]
    sif_meta_data["has matched N"] = sif_meta_data["Paired"].replace(True, "yes").replace(False, "no")
    sif_meta_data["Sex"] = sif_meta_data[SampleAnnotation.sex]
    sif_meta_data["Histology"] = sif_meta_data[SampleAnnotation.histology]

    per_sample_columns = [c for c in meta_data_rows_per_sample if c in meta_data_rows]
    per_patient_columns = [c for c in meta_data_rows if c not in per_sample_columns]
    sif_meta_data = (
        sif_meta_data[[SampleAnnotation.sample, SampleAnnotation.patient] + list(meta_data_rows)]
        .groupby(by=[columns_name])
        .agg({k: list for k in per_sample_columns} | {k: lambda l: np.unique(l)[0] for k in per_patient_columns})
        .reindex(index=columns)
    )
    meta_data = sif_meta_data[meta_data_rows] if meta_data_rows is not None else pd.DataFrame()

    is_all_nan = meta_data.apply(
        lambda row: len([
            x
            for value in row
            for x in (value if isinstance(value, list) else [value])
            if not (isinstance(x, float) and np.isnan(x))
        ]),
        axis=0
    ) == 0
    is_all_nan |= meta_data.isna().all(axis=0)

    meta_data = meta_data[is_all_nan.loc[~is_all_nan].index]

    legend_titles = meta_data.columns
    meta_data_cmaps = {}

    def add_cmap(col, _palette=None, order=None):
        if col not in meta_data.columns:
            return None

        values = []
        for value in meta_data[col].values:
            if isinstance(value, list):
                values += value
            else:
                values.append(value)
        values = [h[0] for h in Counter(values).most_common()]
        if order is not None:
            values = [value for value in order if value in values] + [value for value in values if value not in order]
        if _palette is None:
            cmap = {t: palette.get(t, palette.grey) for t in values}
        else:
            cmap = {value: color for value, color in zip(values, _palette)}
        meta_data_cmaps[col] = cmap

    add_cmap("Sample Type", order=["BM", "EM", "Metastasis", "cfDNA", "T", "P", "N"])
    add_cmap("Material")
    add_cmap("Platform", order=["Agilent", "CCGD", "ICE", "TWIST", "TRACERx", "WGS"])
    add_cmap("has matched N", _palette=[palette.lightgrey, palette.darkgrey], order=["yes", "no"])
    add_cmap("Sex", order=["Female", "Male"])
    add_cmap("Histology")

    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n))
        )
        return new_cmap

    meta_data_cmaps["Contamination"] = truncate_colormap(plt.cm.get_cmap("nipy_spectral"), 0.2, 0.95)
    meta_data_cmaps["Tumor Purity"] = truncate_colormap(plt.cm.get_cmap("nipy_spectral_r"), 0.05, 0.8)

    def meta_data_color(column, value):
        cmap = meta_data_cmaps[column]
        return cmap[value] if isinstance(cmap, dict) else cmap(value)

    ############
    # PLOTTING #
    ############

    # figure layout
    # in units of sample signature exposure square width and height == 1 grid row/column
    pad = 1  # space between subplots
    mut_burden_bar_height = 7
    mut_recurrence_width = 4
    mut_prevalence_width = 3
    cna_prevalence_width = mut_prevalence_width
    total_prevalence_width = 2
    cytoband_wdith = 2
    gene_label_width = 6
    model_annotation_width = 2 if model_annotation is not None else 0
    mut_burden_legend_height = len(mut_burden_legend_cmap)
    mut_legend_height = len(mut_legend_cmap)
    cna_legend_height = len(cna_cmap)
    model_annotation_legend_height = 2
    meta_legend_continuous_height = 4
    meta_legend_height = max(
        [
            len(cmap) if isinstance(cmap, dict) else meta_legend_continuous_height
            for _, cmap in meta_data_cmaps.items()
        ]
    )
    legend_width = 2
    inter_legend_width = 7
    num_legends = max(1, len(legend_titles))
    meta_legend_width = num_legends * (legend_width + inter_legend_width) - inter_legend_width
    inter_heatmap_linewidth = 0.05

    left_of_heatmap_width = (
        mut_recurrence_width
        + pad
        + mut_prevalence_width
        + cna_prevalence_width
        + pad
        + total_prevalence_width
        + pad
        + cytoband_wdith
        + gene_label_width
        + model_annotation_width
        # + 2 * pad
        # + legend_width
    )
    non_heatmap_width = left_of_heatmap_width + 2 * pad + legend_width

    heatmap_height = len(genes)
    heatmap_width = (
        min(len(columns), 10 * max_xfigsize - non_heatmap_width)
        if max_xfigsize is not None
        else len(columns)
    )
    aspect_ratio = heatmap_width / len(columns)
    show_patient_names = (heatmap_width >= len(columns)) and label_columns
    meta_height = len(meta_data.columns) if meta_data is not None else 1

    column_names_height = 5 if show_patient_names else 0

    xsize = left_of_heatmap_width + max(heatmap_width + 2 * pad + legend_width, meta_legend_width)
    ysize = mut_burden_bar_height + pad + heatmap_height + column_names_height
    if meta_data is not None:
        ysize += pad + meta_height + 3 * pad + meta_legend_height

    # xfigsize = xsize / 10 if xfigsize is None else xfigsize
    # yfigsize = ysize / 10 if yfigsize is None else yfigsize
    xfigsize = xsize / 10
    yfigsize = ysize / 10

    ######################
    # PLOTTING FUNCTIONS #
    ######################

    def plot_recurrence(ax, pad=0.1):
        max_comut = max([df.max() for df in mut_protein_data.values() if df is not None] + [1])
        max_xlim = 1.1 * max_comut
        annotations = []
        objects = []

        for row, gene_name in enumerate(genes):
            # paint background rectangle showing area of each gene
            rect = patches.Rectangle(
                xy=(0, row + pad / 2),
                width=max_xlim,
                height=1 - pad,
                facecolor=palette.backgroundgrey,
                edgecolor=None,
                zorder=0.05,
            )
            ax.add_patch(rect)

        for row, (gene_name, df) in enumerate(mut_protein_data.items()):
            if df is None:
                continue

            gene_prot_aggregate = df.sort_values(ascending=True)
            positions = np.linspace(start=row + pad / 2, stop=row + 1 - pad / 2, num=df.size + 2)[1:-1]
            potential_annotations = []
            for pos, ((effect, protein_change), value) in zip(
                positions, gene_prot_aggregate.items()
            ):
                color = palette.get(effect, palette.grey)
                objects.append(
                    ax.hlines(y=pos, xmin=0, xmax=value, color=color, linewidth=0.5)
                )
                objects.append(ax.scatter(value, pos, color=color, s=0.5))
                threshold = snv_comut_threshold if snv_comut_threshold is not None else 5
                if value > threshold + 1 / 6 * (max_comut - threshold):
                    potential_annotations.append(
                        ax.text(
                            value,
                            pos,
                            "  " + protein_change.replace("p.", "") + "  ",
                            fontsize=2,
                            horizontalalignment="right",
                            verticalalignment="top",
                        )
                    )
            if len(potential_annotations) > 1:
                annotations += potential_annotations

        min_xlim = -0.02 * max_xlim  # add enough space for the vertical line at 0 to appear
        ax.set_xlim([min_xlim, max_xlim])
        ax.set_xticks([int(t) for t in ax.get_xticks() if t >= 0])
        ax.set_xlim([min_xlim, max_xlim])
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("top")
        ax.set_xlabel("Recurrence", fontdict=dict(fontsize=5), loc="right")
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=int(np.min([5, max_xlim]))))
        ax.xaxis.set_minor_locator(ticker.FixedLocator([t for t in ax.xaxis.get_minorticklocs() if t >= 0]))
        ax.xaxis.set_minor_formatter(ticker.NullFormatter())

        ax.set_yticks([])
        ax.set_ylim([0, len(genes)])
        ax.yaxis.set_major_locator(ticker.NullLocator())

        ax.tick_params(axis="both", labelsize=4)
        ax.axvline(0, lw=0.7, color=palette.black)
        plotter.grid(ax=ax, axis="x", which="major", zorder=0.5, linewidth=0.5)
        plotter.grid(ax=ax, axis="x", which="minor", zorder=0.1, linewidth=0.3, color=palette.white)
        plotter.no_spines(ax=ax)
        ax.invert_xaxis()

        # ax.set_title("Recurrence", fontsize=8)

        if len(annotations) > 1:
            adjust_text(
                annotations,
                add_objects=objects,
                ax=ax,
                arrowprops=dict(arrowstyle="-", color=palette.black, lw=0.3),
                force_objects=1,
                ha="right",
                va="top",
            )

    def plot_prevalence(ax, pad=0.1):
        max_comut = mut_prevalence_counter.apply(
            lambda c: max(c.values()) if isinstance(c, Counter) else 1
        ).max()
        max_xlim = 1.1 * max_comut

        num_effects = len(effects)

        for row, gene_name in enumerate(genes):
            # paint background rectangle showing area of each gene
            rect = patches.Rectangle(
                xy=(-max_xlim, row + pad / 2),
                width=2 * max_xlim,
                height=1 - pad,
                facecolor=palette.backgroundgrey,
                edgecolor=None,
                zorder=0.05,
            )
            ax.add_patch(rect)

        for row, (gene_name, counter) in enumerate(mut_prevalence_counter.items()):
            if not isinstance(counter, Counter):
                continue

            for i, (effect, color) in enumerate(reversed(mut_cmap.items())):
                if effect in counter:
                    value = counter.get(effect)
                    rect = patches.Rectangle(
                        xy=(0, row + pad / 2 + i * (1 - pad) / num_effects),
                        width=value,
                        height=(1 - pad) / num_effects,
                        facecolor=color,
                        edgecolor=None,
                        zorder=0.9,
                    )
                    ax.add_patch(rect)

        ax.set_xlim([-max_xlim, max_xlim])
        ax.set_xticks([t for t in ax.get_xticks() if t >= 0])
        ax.set_xlim([-max_xlim, max_xlim])
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("top")
        ax.set_xlabel("  SNV", fontdict=dict(fontsize=5), loc="left")
        if max_xlim >= 2:
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=int(min(4, max_xlim))))
            ax.xaxis.set_minor_locator(ticker.FixedLocator([t for t in ax.xaxis.get_minorticklocs() if t >= 0]))
            ax.xaxis.set_minor_formatter(ticker.NullFormatter())

        ax.set_yticks([])
        ax.set_ylim([0, len(genes)])
        ax.yaxis.set_major_locator(ticker.NullLocator())

        ax.tick_params(axis="both", labelsize=4)
        ax.axvline(0, lw=0.7, color=palette.black)
        plotter.grid(ax=ax, axis="x", which="major", zorder=0.5, linewidth=0.5)
        plotter.grid(ax=ax, axis="x", which="minor", zorder=0.1, linewidth=0.3, color=palette.white)
        plotter.no_spines(ax=ax)
        ax.invert_xaxis()

        ax.set_title("Prevalence", fontsize=8)

        cna_ax = ax.twiny()
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("top")

        # get maximum number of amp or del
        max_num_amp = (cna > 0).sum(axis=1).max()
        max_num_del = (cna < 0).sum(axis=1).max()
        max_comut = max(max_num_amp, max_num_del, 1)
        max_xlim = 1.1 * max_comut

        for row, (gene_name, counter) in enumerate(cna_counter.items()):
            non_nan_counter = {k: v for k, v in counter.items() if not np.isnan(k)}
            if not len(non_nan_counter):
                continue

            pos = 0
            for threshold in [high_amp_threshold, low_amp_threshold]:
                if threshold in non_nan_counter:
                    value = non_nan_counter[threshold]
                    rect = patches.Rectangle(
                        xy=(pos, row + 1 / 2),
                        width=value,
                        height=(1 - pad) / 2,
                        facecolor=cna_cmap[threshold],
                        edgecolor=None,
                        zorder=0.9,
                    )
                    cna_ax.add_patch(rect)
                    pos += value
            pos = 0
            for threshold in [high_del_threshold, low_del_threshold]:
                if threshold in non_nan_counter:
                    value = non_nan_counter[threshold]
                    rect = patches.Rectangle(
                        xy=(pos, row + pad / 2),
                        width=value,
                        height=(1 - pad) / 2,
                        facecolor=cna_cmap[threshold],
                        edgecolor=None,
                        zorder=0.9,
                    )
                    cna_ax.add_patch(rect)
                    pos += value

        cna_ax.set_xlim([-max_xlim, max_xlim])
        cna_ax.set_xticks([t for t in cna_ax.get_xticks() if t > 0])
        cna_ax.set_xlim([-max_xlim, max_xlim])
        cna_ax.xaxis.set_ticks_position("top")
        cna_ax.xaxis.set_label_position("top")
        cna_ax.set_xlabel("CNV  ", fontdict=dict(fontsize=5), loc="right")
        cna_ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=int(min(4, max_xlim))))
        cna_ax.xaxis.set_minor_locator(ticker.FixedLocator([t for t in cna_ax.xaxis.get_minorticklocs() if t >= 0]))
        cna_ax.xaxis.set_minor_formatter(ticker.NullFormatter())

        cna_ax.set_yticks([])
        cna_ax.set_ylim([0, len(genes)])
        cna_ax.yaxis.set_major_locator(ticker.NullLocator())

        cna_ax.tick_params(axis="both", labelsize=4)
        cna_ax.axvline(0, lw=0.7, color=palette.black)
        plotter.grid(ax=cna_ax, axis="x", which="major", zorder=0.5, linewidth=0.4)
        plotter.grid(
            ax=cna_ax,
            axis="x",
            which="minor",
            zorder=0.1,
            linewidth=0.3,
            color=palette.white,
        )
        plotter.no_spines(ax=cna_ax)

        def add_percentage_ticks(_ax, xlabel=None):
            percent_ax = _ax.twiny()
            _ax.xaxis.set_ticks_position("top")
            _ax.xaxis.set_label_position("top")
            percent_ax.set_xlim(_ax.get_xlim())
            percent_ax.set_xticks(_ax.get_xticks())
            percent_ax.set_xlim(_ax.get_xlim())
            percent_ax.xaxis.set_ticks_position("bottom")
            percent_ax.xaxis.set_label_position("bottom")
            percent_ax.set_xticklabels(
                [f"{100 * t / len(columns):.2g}%" for t in percent_ax.get_xticks()]
            )
            percent_ax.tick_params(axis="both", labelsize=4)
            percent_ax.xaxis.set_minor_locator(_ax.xaxis.get_minor_locator())
            percent_ax.xaxis.set_minor_formatter(ticker.NullFormatter())
            percent_ax.set_xlabel(xlabel, fontdict=dict(fontsize=6))
            plotter.no_spines(ax=percent_ax)

        _by = "Patients"  # if by == MutA.patient else "Samples"
        add_percentage_ticks(ax, xlabel=f"of {_by}")
        add_percentage_ticks(cna_ax)

        return cna_ax

    def plot_total_prevalence(ax, pad=0.01):
        for row, (gene_name, percentage) in enumerate(total_prevalence_per_gene.items()):
            rect = patches.Rectangle(
                xy=(0, row + pad / 2),
                width=percentage,
                height=1 - pad,
                facecolor=palette.grey,
                edgecolor=None,
                zorder=0.05,
            )
            ax.add_patch(rect)
            ax.text(
                1,
                row + 0.4,
                f"{round(100 * percentage)}%",
                fontsize=3,
                horizontalalignment="right",
                verticalalignment="center",
            )
        ax.set_ylim([0, total_prevalence_per_gene.size])
        ax.set_xlim([0, 1])
        ax.set_xlabel("total", fontdict=dict(fontsize=5), rotation="vertical")
        ax.xaxis.set_label_position("top")
        ax.set_xticks([])
        ax.set_yticks([])
        plotter.no_spines(ax)

    def plot_total_prevalence_overall(ax, pad=0.01):
        rect = patches.Rectangle(
            xy=(0, pad / 2),
            width=total_prevalence_overall,
            height=1 - pad,
            facecolor=palette.grey,
            edgecolor=None,
            zorder=0.05,
        )
        ax.add_patch(rect)
        ax.text(
            1,
            0.4,
            f"{round(100 * total_prevalence_overall)}%",
            fontsize=3,
            horizontalalignment="right",
            verticalalignment="center",
        )
        ax.axhline(y=1, xmin=0, xmax=1, color=palette.black, linewidth=1)
        ax.set_ylim([0, 1])
        ax.set_xlim([0, 1])
        ax.set_xticks([])
        ax.set_yticks([])
        plotter.no_spines(ax)

    def plot_cytoband(ax):
        for y, (gene, cytoband) in enumerate(cytobands.items()):
            ax.text(
                0.5,
                y + 0.4,
                cytoband,
                fontsize=3,
                horizontalalignment="center",
                verticalalignment="center",
            )
        ax.set_ylim([0, cytobands.size])
        ax.set_xlabel("Cytoband", fontdict=dict(fontsize=5), rotation="vertical")
        ax.xaxis.set_label_position("top")
        ax.set_xticks([])
        ax.set_yticks([])
        plotter.no_spines(ax)

    def snv_model_annotation(x, y, width=0.5, height=0.5):
        patch = patches.Ellipse(
            xy=(x + 0.5, y + 0.5),
            width=width,
            height=height,
            facecolor=palette.black,
            edgecolor=None,
        )
        return patch

    def cnv_model_annotation(x, y, width=0.75, height=0.75):
        patch = patches.Rectangle(
            xy=(x + (1 - width) / 2, y + (1 - height) / 2),
            width=width,
            height=height,
            facecolor=palette.white,
            edgecolor=palette.black,
            linewidth=0.25,
        )
        return patch

    def plot_model_annotation(ax):
        for y, (gene, annot) in enumerate(model_annotation.iterrows()):
            x = 0.25  # for better spacing between labels and heatmap
            if annot["cnv"]:
                patch = cnv_model_annotation(x, y)
                ax.add_patch(patch)
            if annot["snv"]:
                patch = snv_model_annotation(x, y)
                ax.add_patch(patch)
        ax.set_xlim([0, 2])
        ax.set_ylim([0, model_annotation.shape[0]])
        # ax.set_xlabel("Cytoband", fontdict=dict(fontsize=5), rotation="vertical")
        # ax.xaxis.set_label_position("top")
        ax.set_xticks([])
        ax.set_yticks([])
        plotter.no_spines(ax)

    def plot_mut_burden(ax, ytickpad=0, fontsize=6):
        mut_burden_data.plot.bar(
            stacked=True,
            width=0.9,
            color=palette,
            ax=ax,
            legend=False,
        )
        ax.set_xlim([-0.5, mut_burden_data.shape[0] - 0.5])
        ax.set_xticks([])
        ax.set_xlabel("")
        ax.set_ylabel(ylabel="Mutation\nBurden", fontdict=dict(fontsize=fontsize))
        ax.yaxis.set_major_formatter(ticker.EngFormatter(sep=""))
        ax.tick_params(axis="y", labelsize=5)
        for tick in ax.yaxis.get_major_ticks():
            tick.set_pad(ytickpad)
        plotter.no_spines(ax)

    def plot_cna_heatmap(ax):
        for (x, y), val in np.ndenumerate(cna.T):
            color = cna_cmap.get(val, palette.white)
            width = 1 - inter_heatmap_linewidth
            height = 1 - inter_heatmap_linewidth * aspect_ratio
            rect = patches.Rectangle(
                xy=(x + (1 - width) / 2, y + (1 - height) / 2),
                width=width,
                height=height,
                facecolor=color,
                edgecolor=None,
            )
            ax.add_patch(rect)

    def plot_mut_heatmap(ax):
        for (x, y), effects in np.ndenumerate(mut.T):
            if isinstance(effects, float):
                continue

            # effective linewidth for the surrounding
            pad = 0.05
            # white background ellipse
            patch = patches.Ellipse(
                xy=(x + 0.5, y + 0.5),
                width=2 / 3 + pad,
                height=2 / 3 + pad,
                facecolor=palette.white,
                edgecolor=None,
            )
            ax.add_patch(patch)

            if len(effects) == 1:
                color = mut_cmap.get(effects[0], palette.grey)
                patch = patches.Ellipse(
                    xy=(x + 0.5, y + 0.5),
                    width=2 / 3,
                    height=2 / 3,
                    facecolor=color,
                    edgecolor=None,
                )
                ax.add_patch(patch)
            else:
                for i, effect in enumerate(reversed(effects)):
                    color = mut_cmap.get(effect, palette.grey)
                    # divide circle into wedges, plotting them clockwise
                    # add small gap to make wedge separation visible
                    theta_offset = 90
                    theta1 = i * 360 / len(effects) + theta_offset
                    theta2 = (i + 1) * 360 / len(effects) + theta_offset
                    theta_pad = np.min([5, (theta2 - theta1) / 10])
                    theta1 += theta_pad
                    theta2 -= theta_pad
                    patch = patches.Wedge(
                        center=(x + 0.5, y + 0.5),
                        r=1 / 3,
                        theta1=theta1,
                        theta2=theta2,
                        facecolor=color,
                        edgecolor=None,
                    )
                    ax.add_patch(patch)

    def plot_heatmap_layout(ax, labelbottom=True):
        ax.set_xlim([0, cna.shape[1]])
        ax.set_ylim([0, cna.shape[0]])
        if labelbottom:
            ax.set_xticks(np.arange(cna.shape[1]) + 0.5)
            ax.set_xticklabels(cna.columns)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.tick_params(axis="x", labelrotation=90)
        ax.tick_params(
            axis="both",
            which="both",
            labelsize=4,
            bottom=False,
            top=False,
            labelbottom=labelbottom,
            left=False,
            right=False,
            labelleft=False,
            pad=0,
        )
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        plotter.no_spines(ax)

    def plot_gene_names(ax):
        ax.set_yticks(np.arange(cna.shape[0]) + 0.5)
        ax.set_yticklabels(genes)
        ax.tick_params(
            axis="y",
            which="both",
            labelsize=4,
            left=False,
            right=False,
            labelleft=True,
            pad=0,
        )
        ax.set_ylabel(None)
        if ground_truth_genes is not None:
            for color, gene_list in ground_truth_genes.items():
                for g, yticklabel in zip(genes, ax.get_yticklabels()):
                    if g in gene_list:
                        yticklabel.set_color(color)
                        yticklabel.set_fontweight("bold")

    def plot_meta_data(ax, labelbottom=True):
        for (x, y), values in np.ndenumerate(meta_data.to_numpy()):
            if not isinstance(values, list):
                values = [values]
            width = 1 - inter_heatmap_linewidth
            height = 1 - inter_heatmap_linewidth * aspect_ratio
            polygons = decompose_rectangle_into_polygons(num=len(values), pad=0.02)
            for i, (val, xy) in enumerate(zip(reversed(values), polygons)):
                ax.fill(
                    x + width * xy[:, 0] + (1 - width) / 2,
                    y + height * xy[:, 1] + (1 - height) / 2,
                    facecolor=meta_data_color(column=legend_titles[y], value=val),
                    edgecolor=None,
                )

        ax.set_xlim([0, meta_data.shape[0]])
        ax.set_ylim([0, meta_data.shape[1]])
        if labelbottom:
            ax.set_xticks(np.arange(meta_data.shape[0]) + 0.5)
            ax.set_xticklabels(meta_data.index)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.tick_params(axis="x", labelrotation=90)
        ax.tick_params(axis="y", labelrotation=0)
        ax.set_yticks(np.arange(meta_data.shape[1]) + 0.5)
        ax.set_yticklabels(legend_titles)
        ax.tick_params(
            axis="both",
            which="both",
            labelsize=4,
            bottom=False,
            top=False,
            labelbottom=labelbottom,
            left=False,
            right=False,
            labelleft=True,
            pad=0,
        )
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.invert_yaxis()
        plotter.no_spines(ax)

    def plot_legend(ax, cmap, names=None, title=None):
        discrete = isinstance(cmap, dict)
        n = len(cmap) if discrete else 2
        if discrete:
            ax.imshow(
                np.arange(n).reshape(n, 1),
                cmap=colors.ListedColormap(list(cmap.values())),
                interpolation="nearest",
                aspect="auto",
            )
        else:
            plt.colorbar(
                mappable=plt.cm.ScalarMappable(cmap=cmap),
                cax=ax
            )
            ax.set_xlim([0, 1])
            ax.set_ylim([0, 1])
        ax.set_xticks([])
        ax.xaxis.set_major_locator(ticker.NullLocator())
        ax.set_yticks(np.arange(n))
        if discrete:
            ax.set_yticklabels(names if names is not None else cmap.keys())
        ax.tick_params(axis="both", labelsize=4, width=0.5)
        ax.yaxis.set_label_position("right")
        ax.yaxis.set_ticks_position("right")
        if title is not None:
            ax.set_title(title, fontsize=4, fontdict=dict(horizontalalignment="left"))
        plotter.set_spines(ax=ax, linewidth=0.5)

    def plot_model_annotation_legend(ax, names=("SNV", "CNV"), title=None):
        # width = 0.5
        # height = 0.5
        # sep = 0.1 / np.sqrt(2)
        # patch = snv_model_annotation(0, 1, width, height, sep)
        patch = snv_model_annotation(0, 1)
        ax.add_patch(patch)
        # patch = cnv_model_annotation(0, 0, width, height, sep)
        patch = cnv_model_annotation(0, 0)
        ax.add_patch(patch)
        ax.set_xlim([-0.5, 1.5])
        ax.set_ylim([0, 2])
        ax.set_xticks([])
        ax.xaxis.set_major_locator(ticker.NullLocator())
        ax.set_yticks(np.arange(2) + 0.5)
        ax.set_yticklabels(reversed(names))
        ax.tick_params(axis="both", labelsize=4, width=0.5)
        ax.yaxis.set_label_position("right")
        ax.yaxis.set_ticks_position("right")
        if title is not None:
            ax.set_title(title, fontsize=4, fontdict=dict(horizontalalignment="left"))
        plotter.no_spines(ax)

    fig = plt.figure(figsize=(xfigsize, yfigsize))
    gs = fig.add_gridspec(ncols=xsize, nrows=ysize)

    x0 = 0  # coordinate origin is at the left of the plot, counting right
    x_recurrence_start = x0  # + pad
    x_recurrence_end = x_recurrence_start + mut_recurrence_width
    x_prevalence_start = x_recurrence_end + pad
    x_prevalence_end = x_prevalence_start + mut_prevalence_width + cna_prevalence_width
    x_total_prevalence_start = x_prevalence_end + pad
    x_total_prevalence_end = x_total_prevalence_start + total_prevalence_width
    x_cytoband_start = x_total_prevalence_end + pad
    x_cytoband_end = x_cytoband_start + cytoband_wdith
    x_model_annotation_start = x_cytoband_end + gene_label_width
    x_model_annotation_end = x_model_annotation_start + model_annotation_width
    x_heatmap_start = x_model_annotation_end
    x_heatmap_end = x_heatmap_start + heatmap_width
    x_meta_start = x_heatmap_start
    x_meta_end = x_heatmap_end
    x_legend_start = x_heatmap_end + 2 * pad
    x_legend_end = x_legend_start + legend_width
    x_meta_legend_start = x_meta_start

    y0 = 0  # coordinate origin is at the top of the plot, counting downwards
    y_burden_start = y0  # + pad
    y_burden_end = y_burden_start + mut_burden_bar_height
    y_burden_legend_start = y_burden_start
    y_burden_legend_end = y_burden_legend_start + mut_burden_legend_height
    y_heatmap_start = y_burden_end + pad
    y_heatmap_end = y_heatmap_start + heatmap_height
    y_total_prevalence_overall_start = y_heatmap_end
    y_total_prevalence_overall_end = y_total_prevalence_overall_start + 1
    y_meta_start = y_heatmap_end + pad
    y_meta_end = y_meta_start + meta_height
    y_mut_legend_start = y_heatmap_start
    y_mut_legend_end = y_mut_legend_start + mut_legend_height
    y_cna_legend_start = y_mut_legend_end + 3* pad
    y_cna_legend_end = y_cna_legend_start + cna_legend_height
    y_model_annotation_legend_start = y_cna_legend_end + 3* pad
    y_model_annotation_legend_end = y_model_annotation_legend_start + model_annotation_legend_height
    y_meta_legend_start = y_meta_end + column_names_height + 3 * pad

    if "mutation burden" in parts_to_plot:
        ax_burden = fig.add_subplot(gs[y_burden_start:y_burden_end, x_heatmap_start:x_heatmap_end])
        plot_mut_burden(ax=ax_burden)
    if "mutation burden legend" in parts_to_plot:
        ax_mut_burden_legend = fig.add_subplot(gs[y_burden_legend_start:y_burden_legend_end, x_legend_start:x_legend_end])
        plot_legend(ax=ax_mut_burden_legend, cmap=mut_burden_legend_cmap)
    if "recurrence" in parts_to_plot:
        ax_recurrence = fig.add_subplot(gs[y_heatmap_start:y_heatmap_end, x_recurrence_start:x_recurrence_end])
        plot_recurrence(ax=ax_recurrence)
    if "prevalence" in parts_to_plot:
        ax_prevalence = fig.add_subplot(gs[y_heatmap_start:y_heatmap_end, x_prevalence_start:x_prevalence_end])
        plot_prevalence(ax=ax_prevalence)
    if "total prevalence" in parts_to_plot:
        ax_total_prevalence = fig.add_subplot(gs[y_heatmap_start:y_heatmap_end, x_total_prevalence_start:x_total_prevalence_end])
        plot_total_prevalence(ax=ax_total_prevalence)
        ax_total_prevalence_overall = fig.add_subplot(gs[y_total_prevalence_overall_start:y_total_prevalence_overall_end, x_total_prevalence_start:x_total_prevalence_end])
        plot_total_prevalence_overall(ax=ax_total_prevalence_overall)
    if "cytoband" in parts_to_plot:
        ax_cytoband = fig.add_subplot(gs[y_heatmap_start:y_heatmap_end, x_cytoband_start:x_cytoband_end])
        plot_cytoband(ax=ax_cytoband)
    if "comutations" in parts_to_plot:
        ax_heatmap = fig.add_subplot(gs[y_heatmap_start:y_heatmap_end, x_heatmap_start:x_heatmap_end])
        plot_cna_heatmap(ax=ax_heatmap)
        plot_mut_heatmap(ax=ax_heatmap)
        plot_heatmap_layout(ax=ax_heatmap, labelbottom=show_patient_names and sif is None)
        if model_annotation_width:
            ax_model_annotation = fig.add_subplot(gs[y_heatmap_start:y_heatmap_end, x_model_annotation_start:x_model_annotation_end])
            plot_model_annotation(ax=ax_model_annotation)
            plot_gene_names(ax=ax_model_annotation)
        else:
            plot_gene_names(ax=ax_heatmap)
    if "comutations legend" in parts_to_plot and mut_legend_height:
        ax_mut_legend = fig.add_subplot(gs[y_mut_legend_start:y_mut_legend_end, x_legend_start:x_legend_end])
        plot_legend(ax=ax_mut_legend, cmap=mut_legend_cmap, title="Short Nucleotide Alteration")
    if "comutations legend" in parts_to_plot and cna_legend_height:
        ax_cna_legend = fig.add_subplot(gs[y_cna_legend_start:y_cna_legend_end, x_legend_start:x_legend_end])
        plot_legend(ax=ax_cna_legend, cmap=cna_cmap, names=cna_names, title="Copy Number Alteration")
    if "comutations legend" in parts_to_plot and model_annotation_width:
        ax_model_annotation_legend = fig.add_subplot(gs[y_model_annotation_legend_start:y_model_annotation_legend_end, x_legend_start:x_legend_end])
        plot_model_annotation_legend(ax=ax_model_annotation_legend, names=model_names, title="Significant by Model")
    if "meta data" in parts_to_plot and meta_height:
        ax_meta = fig.add_subplot(gs[y_meta_start:y_meta_end, x_meta_start:x_meta_end])
        plot_meta_data(ax=ax_meta, labelbottom=show_patient_names)
    if "meta data legend" in parts_to_plot and meta_legend_height:
        xstart = x_meta_legend_start
        for title in legend_titles:
            cmap = meta_data_cmaps[title]
            discrete = isinstance(cmap, dict)
            width = legend_width  # if discrete else 1
            height = len(cmap) if discrete else meta_legend_continuous_height
            xend = xstart + width
            yend = y_meta_legend_start + height
            ax_legend = fig.add_subplot(gs[y_meta_legend_start:yend, xstart:xend])
            plot_legend(ax=ax_legend, cmap=cmap, title=title)
            xstart = xend + inter_legend_width

    name = (
        f"{file_name_prefix}"
        f"comutation_{maf.name}_snv{snv_comut_threshold}|cnv{cnv_comut_threshold}_"
        f"({high_del_threshold}|{low_del_threshold}|{low_amp_threshold}|{high_amp_threshold})"
        f"{file_name_suffix}.pdf"
    )
    plotter.save_figure(
        fig=fig,
        name=name,
        recursive_folder_list=sub_folders,
        bbox_inches="tight",
    )
    plotter.close_figure(fig=fig)
    # print(name)
    return genes, columns

import pandas as pd
from adjustText import adjust_text
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib import colors, patches, ticker
import numpy as np

from src.functional_effect import sort_functional_effects
from src.plotter import Plotter
from src.math import decompose_rectangle_into_polygons
from src.sample_annotation import SampleAnnotation as SA


class ComutPlotter(Plotter):

    def __init__(self, output: str = "./comut.pdf", extra_palette=None) -> None:
        super().__init__(output=output, extra_palette=extra_palette)

    def plot_recurrence(self, ax, mut_protein_data: dict[str, pd.DataFrame], genes: list[str], snv_recurrence_threshold: int = 5, pad=0.1):
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
                facecolor=self.palette.backgroundgrey,
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
                color = self.palette.get(effect, self.palette.grey)
                objects.append(
                    ax.hlines(y=pos, xmin=0, xmax=value, color=color, linewidth=0.5)
                )
                objects.append(ax.scatter(value, pos, color=color, s=0.5))
                if value > snv_recurrence_threshold + 1 / 6 * (max_comut - snv_recurrence_threshold):
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
        ax.axvline(0, lw=0.7, color=self.palette.black)
        self.grid(ax=ax, axis="x", which="major", zorder=0.5, linewidth=0.5)
        self.grid(ax=ax, axis="x", which="minor", zorder=0.1, linewidth=0.3, color=self.palette.white)
        self.no_spines(ax=ax)
        ax.invert_xaxis()

        # ax.set_title("Recurrence", fontsize=8)

        if len(annotations) > 1:
            adjust_text(
                annotations,
                add_objects=objects,
                ax=ax,
                arrowprops=dict(arrowstyle="-", color=self.palette.black, lw=0.3),
                force_objects=1,
                ha="right",
                va="top",
            )

    def plot_prevalence(self, ax, mut_prevalence_counter, cna, cna_counter, num_effects, genes, columns, cnv_cmap, amp_thresholds, del_thresholds, pad=0.1):
        max_comut = mut_prevalence_counter.apply(
            lambda c: max(c.values()) if isinstance(c, Counter) else 1
        ).max()
        max_xlim = 1.1 * max_comut

        effects = sort_functional_effects(
            {
                effect
                for counter in [c for c in mut_prevalence_counter.values if isinstance(c, Counter)]
                for effect, _ in counter.items()
            },
            ascending=False
        )

        for row, gene_name in enumerate(genes):
            # paint background rectangle showing area of each gene
            rect = patches.Rectangle(
                xy=(-max_xlim, row + pad / 2),
                width=2 * max_xlim,
                height=1 - pad,
                facecolor=self.palette.backgroundgrey,
                edgecolor=None,
                zorder=0.05,
            )
            ax.add_patch(rect)

        for row, (gene_name, counter) in enumerate(mut_prevalence_counter.items()):
            if not isinstance(counter, Counter):
                continue

            for i, effect in enumerate(effects):
                if effect in counter:
                    rect = patches.Rectangle(
                        xy=(0, row + pad / 2 + i * (1 - pad) / num_effects),
                        width=counter.get(effect),
                        height=(1 - pad) / num_effects,
                        facecolor=self.palette.get(effect, self.palette.grey),
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
        ax.axvline(0, lw=0.7, color=self.palette.black)
        self.grid(ax=ax, axis="x", which="major", zorder=0.5, linewidth=0.5)
        self.grid(ax=ax, axis="x", which="minor", zorder=0.1, linewidth=0.3, color=self.palette.white)
        self.no_spines(ax=ax)
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
            for threshold in amp_thresholds:
                if threshold in non_nan_counter:
                    value = non_nan_counter[threshold]
                    rect = patches.Rectangle(
                        xy=(pos, row + 1 / 2),
                        width=value,
                        height=(1 - pad) / 2,
                        facecolor=cnv_cmap[threshold],
                        edgecolor=None,
                        zorder=0.9,
                    )
                    cna_ax.add_patch(rect)
                    pos += value
            pos = 0
            for threshold in del_thresholds:
                if threshold in non_nan_counter:
                    value = non_nan_counter[threshold]
                    rect = patches.Rectangle(
                        xy=(pos, row + pad / 2),
                        width=value,
                        height=(1 - pad) / 2,
                        facecolor=cnv_cmap[threshold],
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
        cna_ax.axvline(0, lw=0.7, color=self.palette.black)
        self.grid(ax=cna_ax, axis="x", which="major", zorder=0.5, linewidth=0.4)
        self.grid(
            ax=cna_ax,
            axis="x",
            which="minor",
            zorder=0.1,
            linewidth=0.3,
            color=self.palette.white,
        )
        self.no_spines(ax=cna_ax)

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
            self.no_spines(ax=percent_ax)

        _by = "Patients"  # if by == MutA.patient else "Samples"
        add_percentage_ticks(ax, xlabel=f"of {_by}")
        add_percentage_ticks(cna_ax)

        return cna_ax

    def plot_total_prevalence(self, ax, total_prevalence_per_gene: pd.Series, pad=0.01):
        for row, (gene_name, percentage) in enumerate(total_prevalence_per_gene.items()):
            rect = patches.Rectangle(
                xy=(0, row + pad / 2),
                width=percentage,
                height=1 - pad,
                facecolor=self.palette.grey,
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
        self.no_spines(ax)

    def plot_total_prevalence_overall(self, ax, total_prevalence_overall: float, pad=0.01):
        rect = patches.Rectangle(
            xy=(0, pad / 2),
            width=total_prevalence_overall,
            height=1 - pad,
            facecolor=self.palette.grey,
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
        ax.axhline(y=1, xmin=0, xmax=1, color=self.palette.black, linewidth=1)
        ax.set_ylim([0, 1])
        ax.set_xlim([0, 1])
        ax.set_xticks([])
        ax.set_yticks([])
        self.no_spines(ax)

    def plot_cytoband(self, ax, cytobands: pd.Series):
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
        self.no_spines(ax)

    def snv_model_annotation(self, x, y, width=0.5, height=0.5):
        patch = patches.Ellipse(
            xy=(x + 0.5, y + 0.5),
            width=width,
            height=height,
            facecolor=self.palette.black,
            edgecolor=None,
        )
        return patch

    def cnv_model_annotation(self, x, y, width=0.75, height=0.75):
        patch = patches.Rectangle(
            xy=(x + (1 - width) / 2, y + (1 - height) / 2),
            width=width,
            height=height,
            facecolor=self.palette.white,
            edgecolor=self.palette.black,
            linewidth=0.25,
        )
        return patch

    def plot_model_annotation(self, ax, model_annotation: pd.Series):
        for y, (gene, annot) in enumerate(model_annotation.iterrows()):
            x = 0.25  # for better spacing between labels and heatmap
            if annot["cnv"]:
                patch = self.cnv_model_annotation(x, y)
                ax.add_patch(patch)
            if annot["snv"]:
                patch = self.snv_model_annotation(x, y)
                ax.add_patch(patch)
        ax.set_xlim([0, 2])
        ax.set_ylim([0, model_annotation.shape[0]])
        # ax.set_xlabel("Cytoband", fontdict=dict(fontsize=5), rotation="vertical")
        # ax.xaxis.set_label_position("top")
        ax.set_xticks([])
        ax.set_yticks([])
        self.no_spines(ax)

    def plot_tmb(self, ax, tmb: pd.DataFrame, tmb_threshold=10, ytickpad=0, fontsize=6, aspect_ratio=1):
        if SA.tmb in tmb.columns:
            ax.bar(
                x=np.arange(tmb.shape[0]),
                height=tmb[SA.tmb],
                yerr=tmb[SA.tmb_error],
                align="center",
                color=[self.palette.darkgrey if t < tmb_threshold else self.palette.darkred for t in tmb[SA.tmb]],
                # color=self.palette.darkgrey,
                alpha=1,
                error_kw=dict(ecolor=self.palette.black, lw=aspect_ratio * 0.5),
            )
            ax.set_ylabel(ylabel="TMB [/Mb]", fontdict=dict(fontsize=fontsize))
        else:  # columns are a list of functional effects
            tmb.plot.bar(
                stacked=True,
                width=0.9,
                color=dict(self.palette),
                ax=ax,
                legend=False,
            )
            ax.set_ylabel(ylabel="Mutation\nBurden", fontdict=dict(fontsize=fontsize))
        ax.set_xlim([-0.5, tmb.shape[0] - 0.5])
        ax.set_ylim([0, ax.get_ylim()[1]])
        ax.set_xticks([])
        ax.set_xlabel("")
        ax.yaxis.set_major_formatter(ticker.EngFormatter(sep=""))
        ax.tick_params(axis="y", labelsize=5)
        for tick in ax.yaxis.get_major_ticks():
            tick.set_pad(ytickpad)
        self.no_spines(ax)
        if SA.tmb in tmb.columns:
            # ax.spines["bottom"].set_visible(True)
            self.grid(ax=ax, axis="y", which="major", zorder=0.5)

    def plot_mutsig(self, ax, mutsig, mutsig_cmap, ytickpad=0, fontsize=6):
        fractions = mutsig / mutsig.sum(axis=1).to_numpy()[:, None]
        fractions = fractions[fractions.columns[::-1]]
        fractions.plot.bar(
            stacked=True,
            width=1,
            color=mutsig_cmap,
            ax=ax,
            legend=False,
        )
        ax.set_xlim([-0.5, fractions.shape[0] - 0.5])
        ax.set_ylim([0, 1])
        ax.set_xticks([])
        ax.set_xlabel("")
        ax.set_yticks([0, 1])
        ax.set_ylabel(
            ylabel="Exposure\nFraction",
            fontdict=dict(fontsize=fontsize),
        )
        ax.tick_params(axis="y", labelsize=5)
        for tick in ax.yaxis.get_major_ticks():
            tick.set_pad(ytickpad)
        self.no_spines(ax)

    def plot_cnv_heatmap(self, ax, cnv, cnv_cmap, inter_heatmap_linewidth, aspect_ratio):
        for (x, y), val in np.ndenumerate(cnv.T):
            color = cnv_cmap.get(val, self.palette.white)
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

    def plot_snv_heatmap(self, ax, snv, snv_cmap):
        for (x, y), effects in np.ndenumerate(snv.T):
            if isinstance(effects, float):
                continue

            # effective linewidth for the surrounding
            pad = 0.05
            # white background ellipse
            patch = patches.Ellipse(
                xy=(x + 0.5, y + 0.5),
                width=2 / 3 + pad,
                height=2 / 3 + pad,
                facecolor=self.palette.white,
                edgecolor=None,
            )
            ax.add_patch(patch)

            if len(effects) == 1:
                color = self.palette.get(effects[0], self.palette.grey)
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
                    color = self.palette.get(effect, self.palette.grey)
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

    def plot_heatmap_layout(self, ax, cna, labelbottom=True):
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
        self.no_spines(ax)

    def plot_gene_names(self, ax, genes, ground_truth_genes=None):
        ground_truth_palette = {gene: color for color, gene_list in ground_truth_genes.items() for gene in gene_list} if ground_truth_genes is not None else {}
        for y, gene in enumerate(genes):
            named_color = ground_truth_palette.get(gene)
            color = self.palette.get(named_color, self.palette.black)
            ax.text(
                0.91,
                y + 0.4,
                gene,
                color=color,
                fontsize=4,
                fontweight="bold" if gene in ground_truth_palette else "normal",
                horizontalalignment="right",
                verticalalignment="center",
            )
        ax.set_ylim([0, len(genes)])
        ax.set_xticks([])
        ax.set_yticks([])
        self.no_spines(ax)

    def plot_meta_data(self, ax, meta_data, meta_data_color, legend_titles, inter_heatmap_linewidth, aspect_ratio, labelbottom=True):
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
        self.no_spines(ax)

    def plot_legend(self, ax, cmap, names=None, title=None):
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
        self.set_spines(ax=ax, linewidth=0.5)

    def plot_model_significance(self, ax):
        pass

    def plot_model_annotation_legend(self, ax, names=("SNV", "CNV"), title=None):
        # width = 0.5
        # height = 0.5
        # sep = 0.1 / np.sqrt(2)
        # patch = snv_model_annotation(0, 1, width, height, sep)
        patch = self.snv_model_annotation(0, 1)
        ax.add_patch(patch)
        # patch = cnv_model_annotation(0, 0, width, height, sep)
        patch = self.cnv_model_annotation(0, 0)
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
        self.no_spines(ax)

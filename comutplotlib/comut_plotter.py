import pandas as pd
from adjustText import adjust_text
import matplotlib.pyplot as plt
from matplotlib import colors, patches, ticker
import numpy as np
import scipy.stats as st

from comutplotlib.plotter import Plotter
from comutplotlib.math import decompose_rectangle_into_polygons
from comutplotlib.palette import Palette
from comutplotlib.sample_annotation import SampleAnnotation as SA


class ComutPlotter(Plotter):

    def __init__(self, output: str = "./comut.pdf", extra_palette=None) -> None:
        super().__init__(output=output, extra_palette=extra_palette)

    def plot_recurrence(self, ax, mut_protein_data: dict[str, pd.DataFrame], cna, cna_counter, genes, columns, cnv_cmap, amp_thresholds, del_thresholds, snv_recurrence_threshold: int = 5, pad=0.1):

        def make_snv_recurrence(_ax):
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
                _ax.add_patch(rect)

            for row, gene_name in enumerate(genes):
                df = mut_protein_data.get(gene_name)
                if df is None:
                    continue

                gene_prot_aggregate = df.sort_values(ascending=True)
                positions = np.linspace(start=row + pad / 2, stop=row + 1 - pad / 2, num=df.size + 2)[1:-1]
                potential_annotations = []
                for pos, ((effect, protein_change), value) in zip(positions, gene_prot_aggregate.items()):
                    color = self.palette.get(effect, self.palette.grey)
                    objects.append(_ax.hlines(y=pos, xmin=0, xmax=value, color=color, linewidth=0.5))
                    objects.append(_ax.scatter(value, pos, color=color, s=0.5))
                    if value >= snv_recurrence_threshold:
                        potential_annotations.append(
                            _ax.text(
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

            # min_xlim = -0.02 * max_xlim  # add enough space for the vertical line at 0 to appear
            min_xlim = -max_xlim
            self.set_integer_ticks(ax=_ax, xlim=[min_xlim, max_xlim], xmin=0, n_major=3, n_minor=4)
            _ax.xaxis.set_ticks_position("top")
            _ax.xaxis.set_label_position("top")
            # _ax.set_xlabel("Recurrence", fontdict=dict(fontsize=5), loc="right")

            _ax.set_yticks([])
            _ax.set_ylim([0, len(genes)])
            _ax.yaxis.set_major_locator(ticker.NullLocator())

            _ax.tick_params(axis="both", labelsize=4, pad=1)
            _ax.axvline(0, lw=0.7, color=self.palette.black)
            self.grid(ax=_ax, axis="x", which="major", zorder=0.5, linewidth=0.5)
            self.grid(ax=_ax, axis="x", which="minor", zorder=0.1, linewidth=0.3, color=self.palette.white)
            self.no_spines(ax=_ax)
            _ax.invert_xaxis()

            # ax.set_title("Recurrence", fontsize=8)

            if len(annotations) > 1:
                adjust_text(
                    annotations,
                    # objects=objects,
                    ax=_ax,
                    avoid_self=True,
                    ensure_inside_axes=False,
                    arrowprops=dict(arrowstyle="-", color=self.palette.black, lw=0.3),
                )

        def make_cnv_recurrence(_ax):
            # get maximum number of amp or del
            max_num_amp = (cna > np.min(amp_thresholds)).sum(axis=1).max()
            max_num_del = (cna < np.max(del_thresholds)).sum(axis=1).max()
            max_comut = max(max_num_amp, max_num_del, 1)
            max_xlim = 1.1 * max_comut

            for row, (gene_name, counter) in enumerate(cna_counter.items()):
                non_nan_counter = {k: v for k, v in counter.items() if not np.isnan(k)}
                if not len(non_nan_counter):
                    continue

                present_amp_thresholds = sorted([t for t in amp_thresholds if t in non_nan_counter], reverse=True)
                present_del_thresholds = sorted([t for t in del_thresholds if t in non_nan_counter])

                def make_bars_with_recurrence_text(thresholds, ypad):
                    pos = 0
                    prevalence = {}
                    for threshold in thresholds:
                        if threshold in non_nan_counter:
                            value = non_nan_counter[threshold]
                            rect = patches.Rectangle(
                                xy=(pos, row + pad / 2 + ypad),
                                width=value,
                                height=(1 - pad) / 2,
                                facecolor=cnv_cmap[threshold],
                                edgecolor=None,
                                zorder=0.9,
                            )
                            _ax.add_patch(rect)
                            pos += value
                            prevalence[threshold] = value
                        # Add percentage text for the first and last threshold
                        added_text = False
                        if threshold == amp_thresholds[0] or threshold == del_thresholds[0]:
                            added_text = True
                            _ax.text(
                                0.11 * max_xlim,
                                row + 1 / 4 + pad / 8 + ypad,
                                f"{100 * prevalence.get(threshold, 0) / len(columns):.1f}%",
                                fontsize=1.5,
                                horizontalalignment="left",
                                verticalalignment="center",
                            )
                        if threshold == amp_thresholds[-1] or threshold == del_thresholds[-1]:
                            added_text = True
                            _ax.text(
                                0.91 * max_xlim,
                                row + 1 / 4 + pad / 8 + ypad,
                                f"(+{100 * prevalence.get(threshold, 0) / len(columns):.1f}%)",
                                fontsize=1.5,
                                color=self.palette.grey,
                                horizontalalignment="right",
                                verticalalignment="center",
                            )
                        if not added_text and (threshold == amp_thresholds[1] or threshold == del_thresholds[1]):
                            _ax.text(
                                0.11 * max_xlim,
                                row + 1 / 4 + pad / 8 + ypad,
                                f"{100 * prevalence.get(threshold, 0) / len(columns):.1f}%",
                                fontsize=1.5,
                                horizontalalignment="left",
                                verticalalignment="center",
                            )

                make_bars_with_recurrence_text(present_amp_thresholds, ypad=(1 - pad) / 2)
                make_bars_with_recurrence_text(present_del_thresholds, ypad=0)

            self.set_integer_ticks(ax=_ax, xlim=[-max_xlim, max_xlim], xmin=0.1, n_major=3, n_minor=4)
            _ax.xaxis.set_ticks_position("top")
            _ax.xaxis.set_label_position("top")

            _ax.set_yticks([])
            _ax.set_ylim([0, len(genes)])
            _ax.yaxis.set_major_locator(ticker.NullLocator())

            _ax.tick_params(axis="both", labelsize=4, pad=1)
            _ax.axvline(0, lw=0.7, color=self.palette.black)
            self.grid(ax=_ax, axis="x", which="major", zorder=0.5, linewidth=0.4)
            self.grid(
                ax=_ax,
                axis="x",
                which="minor",
                zorder=0.1,
                linewidth=0.3,
                color=self.palette.white,
            )
            self.no_spines(ax=_ax)

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
                [f"{100 * t / len(columns):.1f}%" for t in percent_ax.get_xticks()]
            )
            percent_ax.tick_params(axis="both", labelsize=4, pad=2)
            percent_ax.xaxis.set_minor_locator(_ax.xaxis.get_minor_locator())
            percent_ax.xaxis.set_minor_formatter(ticker.NullFormatter())
            percent_ax.set_xlabel(xlabel, fontdict=dict(fontsize=4), labelpad=2)
            self.no_spines(ax=percent_ax)

        make_snv_recurrence(_ax=ax)
        cna_ax = ax.twiny()
        label_ax = ax.twiny()
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("top")
        make_cnv_recurrence(_ax=cna_ax)

        _by = "Patients"  # if by == MutA.patient else "Samples"
        add_percentage_ticks(ax, xlabel=f"Fraction of {_by}")
        add_percentage_ticks(cna_ax)

        label_ax.set_xticks([])
        label_ax.set_yticks([])
        label_ax.xaxis.set_ticks_position("top")
        label_ax.xaxis.set_label_position("top")
        self.no_spines(ax=label_ax)

        ax.set_title("Recurrence", fontsize=7)
        ax.set_xlabel("  SNV", fontdict=dict(fontsize=5), loc="left", labelpad=9)
        cna_ax.set_xlabel("CNV  ", fontdict=dict(fontsize=5), loc="right", labelpad=9)
        label_ax.set_xlabel("Number of Patients", fontdict=dict(fontsize=4), labelpad=11)

        return cna_ax

    def plot_total_recurrence(self, ax, total_recurrence_per_gene: pd.DataFrame, pad=0.01, invert_x=False, set_joint_title=None):
        for row, (gene_name, percentage) in enumerate(total_recurrence_per_gene.iterrows()):
            width = percentage["low"]
            x = 1 - width if invert_x else 0
            rect = patches.Rectangle(
                xy=(x, row + pad / 2),
                width=width,
                height=1 - pad,
                facecolor=self.palette.offgrey,
                edgecolor=None,
                zorder=0.05,
            )
            ax.add_patch(rect)
            ax.text(
                -0.08 if invert_x else 1.08,
                row + 0.4,
                f"{round(100 * width)}%",
                fontsize=1.5,
                horizontalalignment="right" if invert_x else "left",
                verticalalignment="center",
                color=self.palette.grey
            )

            width = percentage["high"]
            x = 1 - width if invert_x else 0
            rect = patches.Rectangle(
                xy=(x, row + pad / 2),
                width=width,
                height=1 - pad,
                facecolor=self.palette.grey,
                edgecolor=None,
                zorder=0.05,
            )
            ax.add_patch(rect)
            ax.text(
                0 if invert_x else 1,
                row + 0.4,
                f"{round(100 * width)}%",
                fontsize=3,
                horizontalalignment="left" if invert_x else "right",
                verticalalignment="center",
            )
        ax.set_ylim([0, total_recurrence_per_gene.shape[0]])
        ax.set_xlim([0, 1])
        if set_joint_title is not None:
            if set_joint_title:
                ax.set_xlabel("total", fontdict=dict(fontsize=5), rotation="vertical", x=1.03 if invert_x else -0.03)
        else:
            ax.set_xlabel("total", fontdict=dict(fontsize=5), rotation="vertical")

        ax.xaxis.set_label_position("top")
        ax.set_xticks([])
        ax.set_yticks([])
        self.no_spines(ax)

    def plot_total_recurrence_overall(self, ax, total_recurrence_overall: tuple[dict[str, float], int], pad=0.01, invert_x=False, case_control_spacing=None):
        total_recurrence = total_recurrence_overall[0]
        total = total_recurrence_overall[1]
        percentage = total_recurrence["low"] / total
        x = 1 - percentage if invert_x else 0
        rect = patches.Rectangle(
            xy=(x, pad / 2),
            width=percentage,
            height=1 - pad,
            facecolor=self.palette.offgrey,
            edgecolor=None,
            zorder=0.05,
        )
        ax.add_patch(rect)
        ax.text(
            -0.08 if invert_x else 1.08,
            0.4,
            f"{round(100 * percentage)}%",
            fontsize=1.5,
            horizontalalignment="right" if invert_x else "left",
            verticalalignment="center",
            fontweight="bold",
            color=self.palette.grey
        )

        percentage = total_recurrence["high"] / total
        x = 1 - percentage if invert_x else 0
        rect = patches.Rectangle(
            xy=(x, pad / 2),
            width=percentage,
            height=1 - pad,
            facecolor=self.palette.grey,
            edgecolor=None,
            zorder=0.05,
        )
        ax.add_patch(rect)
        ax.text(
            0 if invert_x else 1,
            0.4,
            f"{round(100 * percentage)}%",
            fontsize=3,
            horizontalalignment="left" if invert_x else "right",
            verticalalignment="center",
            fontweight="bold"
        )
        ax.axhline(y=1, xmin=0, xmax=1, color=self.palette.black, linewidth=1)
        ax.set_ylim([0, 1])
        ax.set_xlim([0, 1])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(
            f"({total_recurrence['high']}/{total})",
            fontdict=dict(fontsize=3),
            labelpad=1,
            loc="center" if case_control_spacing is None else ("right" if case_control_spacing else "left")
        )
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
        ax.set_xticks([])
        ax.set_yticks([])
        self.no_spines(ax)

    def plot_tmb(self, ax, tmb: pd.DataFrame, tmb_threshold=10, ytickpad=0, fontsize=5, aspect_ratio=1, ymin=5 * 1e-1, ymax=1e4, shared_y_ax=None):
        if SA.tmb in tmb.columns:
            if SA.n_vars in tmb.columns and SA.n_bases in tmb.columns:
                # test if tmb is significantly higher than threshold:
                # Jeffrey's prior
                n_vars = tmb[SA.n_vars].to_numpy()
                n_bases = tmb[SA.n_bases].to_numpy()
                dist = st.beta(a=n_vars + 1/2, b=n_bases - n_vars + 1/2)
                yerr_lower = [np.max([ymin, y]) for y in 1e6 * (dist.mean() - dist.ppf(0.05))]
                yerr_upper = [np.min([y, ymax]) for y in 1e6 * (dist.ppf(0.95) - dist.mean())]
                p = dist.cdf(tmb_threshold * 1e-6)
                color = [self.palette.high_tmb if q < 0.05 else self.palette.grey for q in p]
                bar_kwargs = {
                    "yerr": [yerr_lower, yerr_upper],
                    "error_kw": dict(ecolor=self.palette.black, lw=np.sqrt(aspect_ratio) * 0.5),
                }
            elif SA.tmb_error in tmb.columns:
                # test if tmb is significantly higher than threshold:
                # tmb is given in log scale
                p = st.norm(tmb[SA.tmb], np.sqrt(tmb[SA.tmb_error])).cdf(tmb_threshold)
                color = [self.palette.high_tmb if q < 0.05 else self.palette.grey for q in p]
                bar_kwargs = {
                    "yerr": tmb[SA.tmb_error],
                    "error_kw": dict(ecolor=self.palette.black, lw=np.sqrt(aspect_ratio) * 0.5),
                }
            else:
                color = [self.palette.high_tmb if t >= tmb_threshold else self.palette.grey for t in tmb[SA.tmb]]
                bar_kwargs = {}
            perc_high_tmb = color.count(self.palette.high_tmb) / len(color)
            ax.bar(
                x=np.arange(tmb.shape[0]),
                height=tmb[SA.tmb],
                color=color,
                align="center",
                alpha=1,
                **bar_kwargs
            )
            if shared_y_ax is None:
                ax.set_ylabel(ylabel="TMB\n[/Mb]", fontdict=dict(fontsize=fontsize))
        else:  # columns are a list of functional effects
            tmb.plot.bar(
                stacked=True,
                width=0.9,
                color=dict(self.palette),
                ax=ax,
                legend=False,
            )
            if shared_y_ax is None:
                ax.set_ylabel(ylabel="Mutation\nBurden", fontdict=dict(fontsize=fontsize))
            else:
                ax.set_yticks(shared_y_ax.get_yticks())  # set ticks for grid
        ax.set_xlim([-0.5, tmb.shape[0] - 0.5])
        ax.set_yscale("log")
        ax.set_ylim([ymin, ymax])
        ax.set_xticks([])
        ax.set_xlabel("")
        if shared_y_ax is None:
            ax.yaxis.set_minor_locator(ticker.NullLocator())
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: "{:g}".format(y)))
            ax.tick_params(axis="y", labelsize=5)
            for tick in ax.yaxis.get_major_ticks():
                tick.set_pad(ytickpad)
        else:
            ax.yaxis.set_minor_locator(shared_y_ax.yaxis.get_minor_locator())
            ax.yaxis.set_major_formatter(ticker.NullFormatter())
            ax.tick_params(axis="y", length=0)
        self.no_spines(ax)
        if SA.tmb in tmb.columns:
            # ax.spines["bottom"].set_visible(True)
            self.grid(ax=ax, axis="y", which="major", zorder=0.5)
            ax.axhline(y=tmb_threshold, color=self.palette.high_tmb, linewidth=0.5, ls="--", zorder=1)
            ax.axhline(y=tmb[SA.tmb].median(), color=self.palette.darkgrey, linewidth=0.5, ls="--", zorder=0)
            median_tmb = f"{tmb[SA.tmb].median():.1f}"
            median_tmb = " " + median_tmb if shared_y_ax is None else median_tmb + " "
            ax.text(
                x=ax.get_xlim()[1] if shared_y_ax is None else ax.get_xlim()[0],
                y=tmb[SA.tmb].median(),
                s=median_tmb,
                color=self.palette.darkgrey,
                fontsize=4,
                verticalalignment="center",
                horizontalalignment="left" if shared_y_ax is None else "right",
            )
            # create yticklabel on the right for the TMB threshold and label it with perc_high_tmb
            ax2 = ax.twinx()
            if shared_y_ax is not None:
                ax2.yaxis.set_ticks_position("left")
            self.no_spines(ax2)
            ax2.set_ylim(ax.get_ylim())
            ax2.set_yscale("log")
            ax2.set_yticks([tmb_threshold])
            ax2.yaxis.set_minor_locator(ticker.NullLocator())
            ax2.tick_params(axis="y", pad=0, length=0, labelsize=5, labelcolor=self.palette.high_tmb, labelrotation=90)
            ax2.set_yticklabels([" {:.0f}%".format(perc_high_tmb * 100)], fontdict=dict(verticalalignment="bottom"))
            for tick in ax2.yaxis.get_major_ticks():
                tick.set_pad(ytickpad)

    def plot_mutsig(self, ax, mutsig, mutsig_cmap, ytickpad=0, fontsize=5, add_ylabel=True):
        fractions = mutsig / mutsig.sum(axis=1).to_numpy()[:, None]
        fractions = fractions[fractions.columns[::-1]]
        fractions.plot.bar(
            stacked=True,
            width=1,
            color=[mutsig_cmap[c] for c in fractions.columns],
            ax=ax,
            legend=False,
        )
        ax.set_xlim([-0.5, fractions.shape[0] - 0.5])
        ax.set_ylim([0, 1])
        ax.set_xticks([])
        ax.set_xlabel("")
        if add_ylabel:
            ax.set_yticks([0, 1])
            ax.set_ylabel(
                ylabel="Exposure\nFraction",
                fontdict=dict(fontsize=fontsize),
            )
            ax.tick_params(axis="y", labelsize=5)
            for tick in ax.yaxis.get_major_ticks():
                tick.set_pad(ytickpad)
        else:
            ax.set_yticks([])
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

    def plot_meta_data(self, ax, meta_data, meta_data_color, legend_titles, inter_heatmap_linewidth, aspect_ratio, labelbottom=True, add_ylabel=True):
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
        if add_ylabel:
            ax.tick_params(axis="y", labelrotation=0)
            ax.set_yticks(np.arange(meta_data.shape[1]) + 0.5)
            ax.set_yticklabels(legend_titles)
        else:
            ax.set_yticks([])
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

    def plot_legend(self, ax, cmap, names=None, title=None, title_loc="top"):
        discrete = isinstance(cmap, Palette)
        if discrete:
            _cmap = cmap.drop_undefined()
            n = len(_cmap)
            ax.imshow(
                np.arange(n).reshape(n, 1),
                cmap=colors.ListedColormap(list(_cmap.values())),
                interpolation="nearest",
                aspect="auto",
            )
            ax.set_yticks(np.arange(n))
            ax.set_yticklabels(names if names is not None else _cmap.keys())
        else:
            _cmap, _norm = cmap
            plt.colorbar(
                mappable=plt.cm.ScalarMappable(cmap=_cmap, norm=_norm),
                cax=ax
            )
            ax.set_xlim([0, 1])
            ax.set_ylim([_norm.inverse(0), _norm.inverse(1)])
            ax.set_yticks([_norm.inverse(0), _norm.inverse(1)])
        ax.set_xticks([])
        ax.xaxis.set_major_locator(ticker.NullLocator())
        ax.tick_params(axis="both", labelsize=4, width=0.5)
        ax.yaxis.set_label_position("right")
        ax.yaxis.set_ticks_position("right")
        if title is not None:
            if title_loc == "top":
                ax.set_title(title, fontsize=4, loc="left", horizontalalignment="left", pad=5)
            elif title_loc == "bottom":
                ax.set_title(title, fontsize=4, loc="left", horizontalalignment="left", verticalalignment="top", y=0, pad=-5)
            else:
                pass
        self.set_spines(ax=ax, linewidth=0.5)

    def plot_model_significance(self, ax):
        pass

    def plot_model_annotation_legend(self, ax, names=("SNV", "CNV"), title=None, title_loc="top"):
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
            if title_loc == "top":
                ax.set_title(title, fontsize=4, loc="left", horizontalalignment="left", pad=5)
            elif title_loc == "bottom":
                ax.set_title(title, fontsize=4, loc="left", horizontalalignment="left", verticalalignment="top", y=0, pad=-5)
            else:
                pass
        self.no_spines(ax)

from collections import Counter
from itertools import chain
import numpy as np

from comutplotlib.functional_effect import sort_functional_effects


class SNV(object):

    def __init__(self, maf, by):
        self.maf = maf
        self.by = by
        self.df = self.get_df()

    def get_df(self):
        return (
            (self.maf.drop_patient_duplicates() if self.by == self.maf.patient else self.maf)
            .data
            .groupby(by=[self.by, self.maf.gene_name])
            .agg({self.maf.effect: sort_functional_effects})
            .get(self.maf.effect)
            .unstack(self.by)
        )

    def reindex(self, index=None, columns=None):
        if index is not None:
            self.maf.select({self.maf.gene_name: index}, inplace=True)
            self.df = self.df.reindex(index=index)
        if columns is not None:
            self.maf.select({self.by: columns}, inplace=True)
            self.df = self.df.reindex(columns=columns)

    @property
    def has_snv(self):
        return ~self.df.isna()

    @property
    def effects(self):
        effects = set(chain.from_iterable([l for _, l in np.ndenumerate(self.df) if isinstance(l, list)]))
        return sort_functional_effects(effects)

    @property
    def num_effects(self):
        return len(self.effects)
    
    def get_tmb(self):
        pool_as = {
            self.maf.synonymous: {
                self.maf.effect: [
                    self.maf.utr3,
                    self.maf.utr5,
                    self.maf.igr,
                    self.maf.intron,
                    self.maf.silent
                ]
            },
            self.maf.structural: {
                self.maf.effect: [
                    # self.maf.utr3,
                    # self.maf.utr5,
                    self.maf.flank3,
                    self.maf.flank5,
                    self.maf.de_novo_start_in_frame,
                    self.maf.de_novo_start_out_of_frame,
                    self.maf.frame_shift_del,
                    self.maf.frame_shift_ins,
                    # self.maf.igr,
                    # self.maf.intron,
                    self.maf.in_frame_del,
                    self.maf.in_frame_ins,
                    self.maf.linc_rna,
                    # self.maf.missense,
                    self.maf.nonsense,
                    self.maf.nonstop,
                    self.maf.read_through,
                    self.maf.rna,
                    self.maf.silent,
                    self.maf.splice_site,
                    self.maf.start_codon_del,
                    self.maf.start_codon_ins,
                    self.maf.start_codon_snp,
                    self.maf.stop_codon_del,
                    self.maf.stop_codon_ins,
                    self.maf.translation_start_site
                ]
            },
            self.maf.missense: {
                self.maf.effect: [
                    self.maf.missense,
                    # self.maf.in_frame_del,
                    # self.maf.in_frame_ins,
                ]
            },
        }
        tmb = (
            (self.maf.drop_patient_duplicates() if self.by == self.maf.patient else self.maf)
            .pool_annotations(pool_as=pool_as)
            .stratify_mutation_counts(by=[self.by, self.maf.effect], counts_column_name=None)
            .sparse.to_dense()
            .unstack(level=self.maf.effect)
            .fillna(0)
        )
        return tmb[[c for c in tmb.columns if isinstance(c, str)]]

    def get_recurrence(self, index):
        mut_recurrence = (
            self.maf.drop_patient_duplicates()
            .assign_column(
                name="__residue",
                value=self.maf.data[[self.maf.chromosome, self.maf.start_pos, self.maf.protein_change]].apply(
                    lambda row: (
                        f"{row[self.maf.chromosome]}:{row[self.maf.start_pos]}"
                        if isinstance(row[self.maf.protein_change], float) and np.isnan(row[self.maf.protein_change]) or row[self.maf.protein_change] == ""
                        else row[self.maf.protein_change]
                    ),
                    axis=1
                )
            )
            .data[[self.maf.gene_name, self.maf.effect, "__residue"]]
            .replace(np.nan, "NaN")
            .groupby([self.maf.gene_name, self.maf.effect, "__residue"])
            .agg("size")
        )
        return {
            g: (
                mut_recurrence.loc[g]
                if g in mut_recurrence.index.get_level_values(level=self.maf.gene_name)
                else None
            )
            for g in index
        }

    def get_effect_prevalence(self):
        return (
            self.maf.drop_patient_duplicates()
            .data
            .drop_duplicates(subset=[self.maf.patient, self.maf.gene_name, self.maf.effect])
            .groupby(by=[self.maf.gene_name])
            .agg({self.maf.effect: Counter})
            .get(self.maf.effect)
        )

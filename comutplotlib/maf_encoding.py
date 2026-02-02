import pandas as pd

from comutplotlib.mutation_annotation import MutationAnnotation as MutA
# from comutplotlib.nucleotide_context import NucleotideContext


class MAFEncoding(object):

    oncotator = "Oncotator"
    funcotator = "Funcotator"
    tumorportal = "Tumorportal"
    tcga_mc3 = "mc3"

    available_sources = [oncotator, funcotator, tumorportal, tcga_mc3]

    default = funcotator

    def __init__(self, from_source: str):
        assert from_source in self.available_sources
        self.from_source = from_source

    @property
    def column_names_to_default(self):
        if self.from_source == self.oncotator:
            return {
                "HGNC_Ensembl Gene ID": MutA.gene_id,
                "Start_position": MutA.start_pos,
                "End_position": MutA.end_pos,
            }
        elif self.from_source == self.tumorportal:
            return {
                "patient": MutA.patient,
                "gene": MutA.gene_name,
                "classification": MutA.type,
                "type": MutA.effect,
                "chr": MutA.chromosome,
                "pos": MutA.start_pos,
                "ref_allele": MutA.ref_allele,
                "newbase": MutA.alt_allele,
                "context65": MutA.context,
                "cons46": MutA.conservation,
            }
        elif self.from_source == self.tcga_mc3:
            return {
                "CONTEXT": MutA.context
            }
        else:
            return {}

    @property
    def column_names_from_default(self):
        return {default: enc for enc, default in self.column_names_to_default.items()}

    def to_default(self, data: pd.DataFrame):
        # if self.from_source != self.default:
        self.rename_columns(data=data)
        self.rename_chromosomes(data=data)
        self.rename_effects(data=data)
        # self.stringify_context(data=data)
        self.format_mc3(data=data)

    def rename_columns(self, data: pd.DataFrame) -> None:
        return data.rename(columns=self.column_names_to_default, inplace=True)

    @staticmethod
    def rename_chromosomes(data: pd.DataFrame) -> None:
        data.loc[data[MutA.chromosome] == "23", MutA.chromosome] = "X"
        data.loc[data[MutA.chromosome] == "24", MutA.chromosome] = "Y"

    def rename_types(self, data: pd.DataFrame) -> None:
        if self.from_source == self.funcotator:
            data[MutA.type] = data.get(MutA.type).replace(
                {
                    "ONP": MutA.mnv,
                },
            )

    def rename_effects(self, data: pd.DataFrame) -> pd.DataFrame:
        if self.from_source == self.funcotator:
            data[MutA.effect] = data[MutA.effect].replace(
                {
                    "DE_NOVO_START_OUT_FRAME": MutA.de_novo_start_out_of_frame,
                    "DE_NOVO_START_IN_FRAME": MutA.de_novo_start_in_frame,
                    "START_CODON_SNP": MutA.start_codon_snp,
                    "START_CODON_INS": MutA.start_codon_ins,
                    "START_CODON_DEL": MutA.start_codon_del,
                    "COULD_NOT_DETERMINE": None,
                },
            )
        elif self.from_source == self.tumorportal:
            data[MutA.effect] = data[MutA.effect].replace(
                {
                    "Missense": MutA.missense,
                    "Nonsense": MutA.nonsense,
                    "In_frame_Del": MutA.in_frame_del,
                    "In_frame_Ins": MutA.in_frame_ins,
                    "Splice_site": MutA.splice_site,
                    "Splice_Site_Del": MutA.splice_site,
                    "Splice_Site_SNP": MutA.splice_site,
                },
            )
        return data

    # def stringify_context(self, data: pd.DataFrame) -> None:
    #     if self.from_source != self.tumorportal:
    #         return None
    #     data[MutA.context] = data.get(MutA.context).apply(
    #         lambda c: NucleotideContext.from_index(
    #             index=c - 1, ordering=(1, 0, 2)  # contexts are encoded 1-based
    #         ).string
    #     )

    def format_mc3(self, data: pd.DataFrame) -> None:
        if not self.from_source == self.tcga_mc3:
            return None
        data[MutA.sample] = data[MutA.sample_barcode].apply(lambda s: s[:-12] if isinstance(s, str) else s)
        data[MutA.patient] = data[MutA.sample].apply(lambda s: s[:-4] if isinstance(s, str) else s)

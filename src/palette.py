import colorsys
import matplotlib.colors as mc
import numpy as np
import seaborn as sns

from src.mutation_annotation import MutationAnnotation as MutA


class Palette(object):
    (
        blue,
        orange,
        green,
        red,
        violet,
        brown,
        pink,
        grey,
        yellow,
        cyan,
    ) = sns.color_palette()
    (
        brightblue,
        brightorange,
        brightgreen,
        brightred,
        brightviolet,
        brightbrown,
        brightpink,
        brightgrey,
        brightyellow,
        brightcyan,
    ) = sns.color_palette("bright")
    (
        darkblue,
        darkorange,
        darkgreen,
        darkred,
        darkviolet,
        darkbrown,
        darkpink,
        darkgrey,
        darkyellow,
        darkcyan,
    ) = sns.color_palette("dark")
    (
        lightblue,
        lightorange,
        lightgreen,
        lightred,
        lightviolet,
        lightbrown,
        lightpink,
        lightgrey,
        lightyellow,
        lightcyan,
    ) = sns.color_palette("pastel")
    (
        mutedblue,
        mutedorange,
        mutedgreen,
        mutedred,
        mutedviolet,
        mutedbrown,
        mutedpink,
        mutedgrey,
        mutedyellow,
        mutedcyan,
    ) = sns.color_palette("muted")
    (
        colorblindblue,
        colorblindorange,
        colorblindgreen,
        colorblindred,
        colorblindviolet,
        colorblindbrown,
        colorblindpink,
        colorblindgrey,
        colorblindyellow,
        colorblindcyan,
    ) = sns.color_palette("colorblind")
    backgroundgrey = (0.95, 0.95, 0.95)
    white = (1, 1, 1)
    black = (0, 0, 0)

    # for everyone their favo(u)rite spelling:
    gray = grey
    lightgray = lightgrey
    darkgray = darkgrey
    backgroundgray = backgroundgrey

    a, A, c, C, t, T, g, G = sns.color_palette(palette="Paired", n_colors=8)

    def __init__(self, dict: dict = None) -> None:
        self.dict = dict if dict is not None else {
            True: self.lightgray,
            False: self.darkgray,
            np.nan: self.backgroundgray,

            # FUNCTIONAL EFFECTS
            MutA.synonymous: self.darkgray,
            MutA.silent: self.darkgray,
            MutA.utr3: self.gray,
            MutA.utr5: self.gray,
            MutA.flank3: self.gray,
            MutA.flank5: self.gray,
            MutA.intron: self.gray,
            MutA.igr: self.lightgray,
            MutA.missense: self.green,
            MutA.nonsynonymous: self.orange,
            MutA.structural: self.orange,
            MutA.nonsense: self.orange,
            MutA.stop_codon_ins: self.orange,
            MutA.nonstop: self.brown,
            MutA.insertion: self.pink,
            MutA.deletion: self.cyan,
            MutA.frame_shift_ins: self.lightpink,
            MutA.frame_shift_del: self.lightcyan,
            MutA.in_frame_ins: self.pink,
            MutA.in_frame_del: self.cyan,
            MutA.de_novo_start_in_frame: self.lightviolet,
            MutA.de_novo_start_out_of_frame: self.lightviolet,
            MutA.start_codon_snp: self.lightviolet,
            MutA.start_codon_del: self.lightviolet,
            MutA.start_codon_ins: self.lightviolet,
            MutA.splice_site: self.violet,
            MutA.translation_start_site: self.violet,
            MutA.read_through: self.violet,
            MutA.stop_codon_del: self.brown,
            MutA.rna: self.yellow,
            MutA.linc_rna: self.yellow,

            # SELECTION
            "negative selection": self.darkblue,
            "positive selection": self.darkred,
            "case": self.darkpink,
            "control": self.darkcyan,
            "prior case": self.pink,
            "prior control": self.cyan,
            "posterior case": self.darkpink,
            "posterior control": self.darkcyan,

            # BASES
            "A": self.A,
            "C": self.C,
            "G": self.G,
            "T": self.T,

            # SBS MUTATION CHANNELS
            "[>]": self.darkgrey,
            "[C>A]": self.cyan,
            "[C>G]": self.black,
            "[C>T]": self.red,
            "[T>A]": self.grey,
            "[T>C]": self.green,
            "[T>G]": self.pink,
            "[*>*]": self.yellow,

            # INDEL MUTATION CHANNELS
            "[>*]": self.orange,
            "[*>]": self.blue,

            # SAMPLE TYPE
            "BM": self.pink,  # brain metastasis
            "cfDNA": self.violet,  # cell-free DNA
            "EM": self.yellow,  # extra-cranial metastasis
            "ECM": self.yellow,  # extra-cranial metastasis
            "Metastasis": self.brown,
            "N": self.black,  # normal
            "Normal": self.black,  # normal
            "P": self.cyan,  # primary
            "Primary": self.cyan,  # primary
            "N/A": self.grey,

            # SAMPLE MATERIAL
            "FFPE": self.colorblindorange,
            "FF": self.colorblindblue,
            "Plasma": self.colorblindviolet,
            "Blood": self.black,
            "CSF": self.colorblindred,
            "Buffycoat": self.colorblindbrown,
            "Urine": self.colorblindyellow,
            "Cell line": self.colorblindgreen,

            # SEX GENOTYPE
            "XX": self.darkpink,
            "Female": self.darkpink,
            "XY": self.darkcyan,
            "Male": self.darkcyan,
            "NA": self.grey,
            "unknown": self.grey,

            # HISTOLOGIES
            "AML": self.darkbrown,  # akute myeloid leukemia
            "CLL": self.lightbrown,  # chronic lymphocytic leukemia
            "LEUK": self.brightbrown,  # leukemia

            "BRCA": self.brightred,  # breast cancer

            "DLBCL": self.red,  # diffuse large B-cell lymphoma
            "CARC": self.lightred,  # carcinoid
            "ECRI": self.darkred,  # endocrine cancer

            "CNC": self.mutedred,  # central neurocytoma
            "ESNB": self.colorblindpink,  # Esthesioneuroblastoma
            "GBM": self.lightpink,  # glioblastoma multiforme
            "GLIO": self.darkpink,  # glioma
            "MED": self.mutedpink,  # medulloblastoma
            "MEN": self.brightpink,  # meningioma
            "NB": self.mutedbrown,  # neuroblastoma
            "PCNSL": self.pink,  # primary central nervous system lymphoma
            "RHAB": self.mutedgrey,  # rhabdoid tumor

            "HNC": self.brightyellow,  # head and neck cancer
            "HNSC": self.brightyellow,  # head and neck cancer
            "LACA": self.lightyellow,  # laryngeal cancer
            "THCA": self.darkyellow,  # thyroid cancer

            "GI": self.brightcyan,  # gastro-intestinal
            "ESO": self.lightcyan,  # esophageal adenocarcinoma
            "CHOL": self.colorblindcyan,  # Cholangiocarcinoma (Bile duct cancer)
            "BLCA": self.cyan,  # bladder cancer
            "CRC": self.darkcyan,  # colorectal cancer

            "KIRC": self.brightgreen,  # kidney cancer
            "RCC": self.brightgreen,  # renal cell carninoma

            "LUAD": self.brightblue,  # lung adenocarcinoma
            "LUCA": self.lightblue,  # lung cancer

            "MEL": self.brightorange,  # melanoma

            "MM": self.darkgray,  # multiple myeloma (bone marrow)

            "SARC": self.gray,  # sarcoma

            "PRAD": self.lightviolet,  # prostrate adenocarcinoma
            "TECA": self.colorblindviolet,  # testicular cancer
            "GynOnc": self.brightviolet,  # ovarian, endometrial, uterine, vaginal cancer
            "OV": self.violet,  # ovarian cancer
            "UCEC": self.darkviolet,  # endometrial cancer

            "hd": self.white,  # healthy donor

            # PLATFORM
            "Agilent": self.lightred,
            "CCGD": self.lightorange,
            "ICE": self.lightcyan,
            "TWIST": self.lightgreen,
            "WES": self.lightorange,
            "WGS": self.lightblue,
            "TRACERx": self.lightviolet,
            "ULP": self.lightgrey,
        }

    def __getitem__(self, item):
        return self.dict.get(item)

    def __setitem__(self, key, value):
        self.dict[key] = value

    def copy(self):
        return Palette(dict=self.dict)

    @staticmethod
    def adjust_lightness(color, amount):
        """from
        https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
        """
        try:
            c = mc.cnames[color]
        except:
            c = color
        c = colorsys.rgb_to_hls(*mc.to_rgb(c))
        return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

    @staticmethod
    def make_diverging_palette(color, n_colors: int) -> list:
        palette1 = sns.light_palette(color, n_colors=np.ceil(n_colors / 2) + 2)[1:-1]
        palette2 = sns.dark_palette(color, n_colors=np.floor(n_colors / 2) + 3)[3:]
        return [p for p in palette1] + [p for p in reversed(palette2)]

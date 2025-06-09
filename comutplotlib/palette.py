from collections import Counter, defaultdict, UserDict
import colorsys
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from comutplotlib.functional_effect import better_effect_legend
from comutplotlib.mutation_annotation import MutationAnnotation as MutA


class Palette(UserDict):
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
    (
        deepblue,
        deeporange,
        deepgreen,
        deepred,
        deepviolet,
        deepbrown,
        deeppink,
        deepgrey,
        deepyellow,
        deepcyan,
    ) = sns.color_palette("deep")

    offgrey = (0.90, 0.90, 0.90)
    backgroundgrey = (0.95, 0.95, 0.95)
    white = (1, 1, 1)
    black = (0, 0, 0)

    grassgreen = (0.0, 160 / 255, 0.0)
    teal = (102 / 255, 194 / 255, 165 / 255)
    salmon = (252 / 255, 141 / 255, 98 / 255)

    # for everyone their favo(u)rite spelling:
    gray = grey
    lightgray = lightgrey
    darkgray = darkgrey
    backgroundgray = backgroundgrey

    a, A, c, C, t, T, g, G = sns.color_palette(palette="Paired", n_colors=8)

    @classmethod
    def normalizeRGB(cls, r: int, g: int, b: int):
        return r / 255, g / 255, b / 255

    @classmethod
    def mix(cls, c1, c2, weight=0.5):
        return tuple(weight * np.array(c1) + (1 - weight) * np.array(c2))

    @classmethod
    def from_hash(cls, hash):
        return Palette({entry.split(":")[0]: [float(c) for c in entry.split(":")[1].split(",")] for entry in hash.split(";")})

    def __init__(self, dict: dict = None) -> None:
        self.high_tmb = self.mix(self.darkred, self.grey)
        super().__init__(
            dict if dict is not None else {
                True: self.lightgrey,
                False: self.darkgrey,
                "yes": self.lightgrey,
                "no": self.darkgrey,
                "NA": self.grey,
                np.nan: self.white,
                "nan": self.white,
                "unknown": self.white,

                # FUNCTIONAL EFFECTS
                MutA.synonymous: self.darkgrey,
                MutA.silent: self.darkgrey,
                MutA.utr3: self.grey,
                MutA.utr5: self.grey,
                MutA.flank3: self.grey,
                MutA.flank5: self.grey,
                MutA.intron: self.grey,
                MutA.igr: self.lightgrey,
                MutA.missense: self.grassgreen,
                MutA.nonsynonymous: self.black,
                MutA.structural: self.black,
                MutA.nonsense: self.black,
                MutA.stop_codon_ins: self.black,
                MutA.nonstop: self.darkbrown,
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

                # CASE - CONTROL
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

                # MUTATIONAL SIGNATURES
                "Ageing": self.brown,
                "Smoking": self.cyan,
                "UV": self.orange,
                "Chemotherapy": self.green,

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
                "FFPE": self.colorblindyellow,
                "FF": self.colorblindcyan,
                "Plasma": self.colorblindviolet,
                "Blood": self.colorblindgrey,
                "CSF": self.colorblindpink,
                "Buffycoat": self.colorblindbrown,
                "Urine": self.colorblindorange,
                "Cell line": self.colorblindred,

                # SEX GENOTYPE
                "XX": self.darkpink,
                "Female": self.darkpink,
                "XY": self.darkcyan,
                "Male": self.darkcyan,

                # HISTOLOGIES
                "AML": self.darkbrown,  # akute myeloid leukemia
                "Acute myeloid leukemia": self.darkbrown,
                "CLL": self.lightbrown,  # chronic lymphocytic leukemia
                "Chromic lymphocytic leukemia": self.lightbrown,
                "LEUK": self.brightbrown,  # leukemia
                "Leukemia": self.brightbrown,

                "BRCA": self.brightred,  # breast cancer
                "Breast cancer": self.brightred,

                "DLBCL": self.red,  # diffuse large B-cell lymphoma
                "Diffuse large B-cell lymphoma": self.red,
                "CARC": self.lightred,  # carcinoid
                "Carcinoid": self.lightred,
                "ECRI": self.darkred,  # endocrine cancer
                "Endocrine cancer": self.darkred,

                "CNC": self.mutedred,  # central neurocytoma
                "Central neurocytoma": self.mutedred,
                "ESNB": self.colorblindpink,  # Esthesioneuroblastoma
                "Esthesioneuroblastoma": self.colorblindpink,
                "GBM": self.lightpink,  # glioblastoma multiforme
                "Glioblastoma multiforme": self.lightpink,
                "GLIO": self.darkpink,  # glioma
                "Glioma": self.darkpink,
                "MED": self.mutedpink,  # medulloblastoma
                "Medulloblastoma": self.mutedpink,
                "MEN": self.brightpink,  # meningioma
                "Meningioma": self.brightpink,
                "NB": self.mutedbrown,  # neuroblastoma
                "Neuroblastoma": self.mutedbrown,
                "PCNSL": self.pink,  # primary central nervous system lymphoma
                "Primary central nervous system lymphoma": self.pink,
                "RHAB": self.mutedgrey,  # rhabdoid tumor
                "Rhabdoid tumor": self.mutedgrey,

                "HNC": self.brightyellow,  # head and neck cancer
                "Head and neck cancer": self.brightyellow,
                "HNSC": self.brightyellow,  # head and neck cancer
                "Head and neck squamous cell carcinoma": self.brightyellow,
                "LACA": self.lightyellow,  # laryngeal cancer
                "Laryngeal cancer": self.lightyellow,
                "THCA": self.darkyellow,  # thyroid cancer
                "Thyroid cancer": self.darkyellow,

                "GI": self.brightcyan,  # gastro-intestinal
                "Gastro-intestinal": self.brightcyan,
                "ESO": self.lightcyan,  # esophageal adenocarcinoma
                "Esophageal adenocarcinoma": self.lightcyan,
                "CHOL": self.colorblindcyan,  # Cholangiocarcinoma (Bile duct cancer)
                "Cholangiocarcinoma": self.colorblindcyan,
                "BLCA": self.cyan,  # bladder cancer
                "Bladder cancer": self.cyan,
                "CRC": self.darkcyan,  # colorectal cancer
                "Colorectal cancer": self.darkcyan,

                "KIRC": self.brightgreen,  # renal cell carcinoma (clear cell)
                "KIRP": self.lightgreen,  # renal cell carcinoma (papillary)
                "KICH": self.darkgreen,  # renal cell carcinoma (chromophobe)
                "RCC": self.brightgreen,  # renal cell carcinoma
                "Renal cell carcinoma": self.brightgreen,
                "Kidney cancer": self.brightgreen,

                "LUAD": self.brightblue,  # lung adenocarcinoma
                "Lung adenocarcinoma": self.brightblue,
                "LUCA": self.lightblue,  # lung cancer
                "Lung cancer": self.lightblue,
                "LUSC": self.darkblue,  # lung squamous cell carcinoma
                "Lung squamous cell carcinoma": self.darkblue,

                "MEL": self.brightorange,  # melanoma
                "Melanoma": self.brightorange,

                "MM": self.darkgrey,  # multiple myeloma (bone marrow)
                "Multiple myeloma": self.darkgrey,

                "SARC": self.grey,  # sarcoma
                "Sarcoma": self.grey,

                "PRAD": self.lightviolet,  # prostrate adenocarcinoma
                "Prostate adenocarcinoma": self.lightviolet,
                "TECA": self.colorblindviolet,  # testicular cancer
                "Testicular cancer": self.colorblindviolet,
                "GynOnc": self.brightviolet,  # ovarian, endometrial, uterine, vaginal cancer
                "OV": self.violet,  # ovarian cancer
                "Ovarian cancer": self.violet,
                "UCEC": self.darkviolet,  # endometrial cancer
                "Endometrial cancer": self.darkviolet,

                "hd": self.white,  # healthy donor
                "Healthy donor": self.white,

                # Hormone Receptor Status
                "HR+": self.yellow,
                "HER2+": self.red,
                "HR+/HER2+": self.yellow,
                "HR-/HER2+": self.red,
                "HR+/HER2-": self.blue,
                "HR-/HER2-": self.darkviolet,
                "TN": self.darkviolet,
                "pos": self.normalizeRGB(102, 194, 165),  # teal
                "neg": self.normalizeRGB(252, 141,  98),  # salmon
                "equiv": self.grey,

                # SEQUENCING PLATFORM
                "Agilent": self.lightcyan,
                "Agilent BI": self.lightcyan,
                "Agilent CCGD": self.lightgreen,
                "ICE": self.lightyellow,
                "TWIST": self.lightred,
                "NimGen hg18": self.lightblue,
                "NimGen v2": self.lightyellow,
                "NimGen v3": self.mutedyellow,
                "NimGen VCRome": self.lightviolet,
                "NimGen VCRome-PKv1": self.mutedviolet,
                "Agilent SS 38Mb": self.lightorange,
                "Agilent SS 50Mb": self.mutedorange,
                "PanCanPanel": self.lightyellow,
                "TRACERx": self.lightviolet,

                "WES": self.lightorange,
                "WGS": self.lightblue,
                "Panel": self.lightyellow,
                "ULP": self.lightgrey,
            }
        )

    def drop_undefined(self):
        return Palette({
            v: c for v, c, in self.items()
            if not (isinstance(v, float) and np.isnan(v) or v in ["nan", "unknown"])
        })

    def hash(self):
        return ";".join([f"{k}:" + ",".join([str(c) for c in self[k]]) for k in self.keys()])

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

    def get_snv_cmap(self, snv):
        # Some functional effects are painted with the same color, which messes up
        # the legend, so we have to construct the legend color map separately.
        mut_cmap = {}
        inv_mut_cmap = defaultdict(list)
        for effect in snv.effects:
            color = self.get(effect, self.grey)
            mut_cmap[effect] = color
            inv_mut_cmap[color].append(better_effect_legend(effect))
        snv_cmap = {
            " / ".join(e_list): c for c, e_list in inv_mut_cmap.items()
        }
        return Palette(snv_cmap)

    def get_cnv_cmap(self, cnv, full=False):
        amp_color = self.red
        del_color = self.blue
        _cnv_cmap = {
            cnv.high_amp_threshold: self.adjust_lightness(amp_color, 1.17),
            cnv.mid_amp_threshold: self.adjust_lightness(amp_color, 1.2),
            cnv.low_amp_threshold: self.adjust_lightness(amp_color, 1.78),
            cnv.baseline: self.backgroundgrey,
            cnv.low_del_threshold: self.adjust_lightness(del_color, 2.03),
            cnv.mid_del_threshold: self.adjust_lightness(del_color, 1.2),
            cnv.high_del_threshold: self.adjust_lightness(del_color, 1.07),
        }
        _cnv_names = [
            "Amplification",
            "Amplification",
            "Gain",
            "Baseline",
            "Loss",
            "Deletion",
            "Deletion",
        ]
        cnv_cmap = {}
        cnv_names = []
        for name, (key, color) in zip(_cnv_names, _cnv_cmap.items()):
            if key in np.unique(cnv.df) or full:
                cnv_cmap[key] = color
                cnv_names.append(name)
        return Palette(cnv_cmap), cnv_names

    def get_tmb_cmap(self, tmb):
        tmb_cmap = {
            better_effect_legend(effect): self.get(effect, self.grey)
            for effect in reversed(tmb.columns)
        } if tmb is not None else {}
        return Palette(tmb_cmap)

    def get_mutsig_cmap(self, mutsig):
        # The listed order determines the order in which the signatures are plotted
        signature_sets = {
            "clock-like": ["SBS1", "SBS5"],
            "PolE/D/H": ["SBS9", "SBS10a", "SBS10b", "SBS10c", "SBS10d", "DBS3"],
            "MMR": [
                "SBS3",
                "SBS6",
                "SBS14",
                "SBS15",
                "SBS20",
                "SBS21",
                "SBS26",
                "SBS30",
                "SBS36",
                "SBS44",
                "DBS7",
                "DBS10",
                "DBS13",
                "ID6",
                "ID7",
                "ID8",  # TOP2A
                "ID17",  # TOP2A
            ],
            "Other": [
                "SBS18",  # Oxog
                "SBS22a",  # Aristolochic acid
                "SBS22b",  # Aristolochic acid
                "SBS24",  # aflatoxin
                "SBS42",  # haloalkane
                "SBS85",  # ind eff of AID
                "SBS88",  # colibactin
                "SBS90",  # duocarmycin
                "SBS99",  # melphalan
                "DBS20",  # Aristolochic acid
                "ID23",  # Aristolochic acid
                "ID18",  # e.coli
            ],
            "APOBEC": ["SBS2", "SBS13"],
            "Smoking": ["SBS4", "SBS29", "SBS92", "DBS2", "ID3"],
            "UV": ["SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38", "DBS1", "ID13"],
            "Treatment": [
                "SBS11",  # Temozolomide
                "SBS25",  # Chemotherapy
                "SBS31",  # Platinum chemotherapy
                "SBS32",  # Azathioprine
                "SBS35",  # Platinum chemotherapy
                "SBS86",  # Unknown chemotherapy
                "SBS87",  # Thiopurine chemotherapy
                "DBS5",  # Platinum chemotherapy
            ],
            "Error": [
                "SBS27",
                "SBS43", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49",
                "SBS50", "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57", "SBS58", "SBS59",
                "SBS60",
                "SBS95",
                "DBS14"
            ],
            "Unknown": []
        }
        mutsigset_palette = {
            "clock-like": self.lightbrown,
            "APOBEC": self.lightred,
            "MMR": self.lightgreen,
            "PolE/D/H": self.lightcyan,
            "Other": self.lightviolet,
            "UV": self.lightorange,
            "Smoking": self.lightblue,
            "Treatment": self.lightpink,
            "Error": self.lightgrey
        }
        signature_colors = {color: signature_sets[key] for key, color in mutsigset_palette.items()}
        for color, signatures in signature_colors.items():
            for s, c in zip(signatures, self.make_diverging_palette(color=color, n_colors=len(signatures))):
                mutsigset_palette[s] = c
        mutsig_cmap = {
            sig: mutsigset_palette.get(sig, color)
            for sig, color in zip(mutsig.columns, sns.color_palette("husl", n_colors=len(mutsig.columns)))
        } if mutsig is not None else {}
        return Palette(mutsig_cmap)

    def get_meta_cmaps(self, meta):
        meta_cmaps = {}

        def add_cmap(col, _palette=None, order=None):
            if col not in meta.rows:
                return None

            values = []
            for value in meta.df[col].values:
                if isinstance(value, list):
                    values += value
                else:
                    values.append(value)

            values = [h[0] for h in Counter(values).most_common()]
            if order is not None:
                _values = [v for v in order if v in values]
                if _palette is not None:
                    _palette = _palette[:len(_values)]
                _values += [v for v in values if v not in order]
                values = _values

            if _palette is None:
                cmap = {t: self.get(t, self.grey) for t in values}
            else:
                cmap = {v: c for v, c in zip(values, _palette)} | {v: self.get(v, self.grey) for v in values[len(_palette):]}
            meta_cmaps[col] = Palette(cmap)

        def add_cont_cmap(col, cmap, minval, maxval, n=100, center=None):
            if col not in meta.rows:
                return None

            if center is None:
                norm = mc.Normalize(vmin=minval, vmax=maxval)
                truncated_cmap = mc.LinearSegmentedColormap.from_list(
                    name=f"trunc({cmap.name},{norm(minval):.2f},{norm(maxval):.2f})",
                    colors=cmap(norm(np.linspace(minval, maxval, n)))
                )
            else:
                norm = mc.TwoSlopeNorm(vmin=minval, vcenter=center, vmax=maxval)
                truncated_cmap = mc.LinearSegmentedColormap.from_list(
                    name=f"trunc({cmap.name},{norm(minval):.2f},{norm(maxval):.2f}, center={center})",
                    colors=cmap(norm(np.linspace(minval, maxval, n)))
                )
            meta_cmaps[col] = (truncated_cmap, norm)

        add_cmap("Sample Type", order=["BM", "EM", "Metastasis", "cfDNA", "T", "P", "N"])
        add_cmap("Material", order=["FF", "FFPE"])
        add_cmap("Platform", order=["Agilent BI", "Agilent CCGD", "ICE", "TWIST", "TRACERx", "WES", "WGS"])
        add_cmap("has matched N", order=["yes", "no"])
        add_cmap("has metastasis", order=["yes", "no"])
        add_cmap("Sex", _palette=[self.normalizeRGB(251, 180, 196), self.normalizeRGB(170, 245, 232)], order=["Female", "Male"])
        add_cmap("HR Status", order=["HR+", "HER2+", "HR+/HER2+", "HR-/HER2+", "HR+/HER2-", "HR-/HER2-", "TN", "N/A", "NA", "pos", "neg", "unknown"])
        add_cmap("ER status", order=["pos", "neg", "unknown"])
        add_cmap("PR status", order=["pos", "neg", "unknown"])
        add_cmap("HER2 status", order=["pos", "neg", "unknown"])
        add_cmap("HER2 Status", order=["pos", "neg", "unknown"])
        add_cmap("ESR1 status", order=["pos", "neg", "unknown"])
        add_cmap("PGR status", order=["pos", "neg", "unknown"])
        add_cmap("ERBB2 status", order=["pos", "neg", "unknown"])
        add_cmap("Histology")
        for col in ["Chemo Tx", "XRT", "Targeted Tx", "Hormone Tx", "Immuno Tx ICI", "ADC"]:
            add_cmap(col, _palette=[self.normalizeRGB(105, 200, 219), self.lightgrey, self.backgroundgrey], order=["yes", "no", "unknown"])
        add_cmap("WGD", _palette=self.make_diverging_palette(self.violet, n_colors=5*4)[::5], order=list(range(4)))
        add_cont_cmap("Contamination", plt.cm.get_cmap("BuPu"), 0.0, 0.05)
        add_cont_cmap("Tumor Purity", plt.cm.get_cmap("plasma_r"), 0, 1)
        add_cont_cmap("Ploidy", plt.cm.get_cmap("PiYG"), 1, 6, center=2)
        add_cont_cmap("Subclonal Fraction", plt.cm.get_cmap("RdPu"), 0, 0.5)
        add_cont_cmap("Age at BM Dx", plt.cm.get_cmap("bone_r"), 0, 100)
        add_cont_cmap("Age at P Dx", plt.cm.get_cmap("bone_r"), 0, 100)

        # For all custom columns:
        for col in meta.rows:
            if col not in meta_cmaps.keys():
                add_cmap(col)

        return {k: meta_cmaps[k] for k in meta.rows if k in meta_cmaps.keys()}

    @staticmethod
    def condense(cmaps):
        continuous_cmaps = {k: p for k, p in cmaps.items() if not isinstance(p, Palette)}
        if "Age at BM Dx" in continuous_cmaps and "Age at P Dx" in continuous_cmaps:
            continuous_cmaps["Age at Dx"] = continuous_cmaps["Age at P Dx"]
            del continuous_cmaps["Age at BM Dx"]
            del continuous_cmaps["Age at P Dx"]
        condenseable_cmaps = {k: p for k, p in cmaps.items() if isinstance(p, Palette)}
        grouped_palettes = defaultdict(list)
        for k, p in condenseable_cmaps.items():
            grouped_palettes[p.hash()].append(k)

        condensed_cmaps = {"\n".join(k_list): Palette.from_hash(h) for h, k_list in grouped_palettes.items()}
        return continuous_cmaps | condensed_cmaps

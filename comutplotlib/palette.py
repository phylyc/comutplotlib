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
    backgroundgrey = (0.95, 0.95, 0.95)
    white = (1, 1, 1)
    black = (0, 0, 0)

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
    def from_hash(cls, hash):
        return Palette({entry.split(":")[0]: [float(c) for c in entry.split(":")[1].split(",")] for entry in hash.split("_")})

    def __init__(self, dict: dict = None) -> None:
        super().__init__(
            dict if dict is not None else {
                True: self.lightgray,
                False: self.darkgray,
                "yes": self.lightgray,
                "no": self.darkgray,
                "NA": self.grey,
                np.nan: self.white,
                "nan": self.white,
                "unknown": self.white,

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

                "KIRC": self.brightgreen,  # renal cell carcinoma (clear cell)
                "KIRP": self.lightgreen,  # renal cell carcinoma (papillary)
                "RCC": self.brightgreen,  # renal cell carcinoma

                "LUAD": self.brightblue,  # lung adenocarcinoma
                "LUCA": self.lightblue,  # lung cancer
                "LUSC": self.darkblue,  # lung squamous cell carcinoma

                "MEL": self.brightorange,  # melanoma

                "MM": self.darkgray,  # multiple myeloma (bone marrow)

                "SARC": self.gray,  # sarcoma

                "PRAD": self.lightviolet,  # prostrate adenocarcinoma
                "TECA": self.colorblindviolet,  # testicular cancer
                "GynOnc": self.brightviolet,  # ovarian, endometrial, uterine, vaginal cancer
                "OV": self.violet,  # ovarian cancer
                "UCEC": self.darkviolet,  # endometrial cancer

                "hd": self.white,  # healthy donor

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

                # PLATFORM
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
                "WES": self.lightorange,
                "WGS": self.lightblue,
                "TRACERx": self.lightviolet,
                "ULP": self.lightgrey,
            }
        )

    def drop_undefined(self):
        return Palette({
            v: c for v, c, in self.items()
            if not (isinstance(v, float) and np.isnan(v) or v in ["nan", "unknown"])
        })

    def hash(self):
        # sorted_keys = sorted(self.keys())
        return "_".join([f"{k}:" + ",".join([str(c) for c in self[k]]) for k in self.keys()])

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

    def get_snv_cmap(self, data):
        # Some functional effects are painted with the same color, which messes up
        # the legend, so we have to construct the legend color map separately.
        mut_cmap = {}
        inv_mut_cmap = defaultdict(list)
        for effect in data.snv.effects:
            color = self.get(effect, self.grey)
            mut_cmap[effect] = color
            inv_mut_cmap[color].append(better_effect_legend(effect))
        snv_cmap = {
            " / ".join(e_list): c for c, e_list in inv_mut_cmap.items()
        }
        return Palette(snv_cmap)

    def get_cnv_cmap(self, data):
        amp_color = self.red
        del_color = self.blue
        _cnv_cmap = {
            data.high_amp_threshold: self.adjust_lightness(amp_color, 1.2),
            data.mid_amp_threshold: self.adjust_lightness(amp_color, 1.2),
            data.low_amp_threshold: self.adjust_lightness(amp_color, 1.8),
            data.baseline: self.backgroundgrey,
            data.low_del_threshold: self.adjust_lightness(del_color, 2.05),
            data.mid_del_threshold: self.adjust_lightness(del_color, 1.2),
            data.high_del_threshold: self.adjust_lightness(del_color, 1.1),
        }
        _cnv_names = [
            "High Amplification",
            "Amplification",
            "Low Amplification",
            "Baseline",
            "Shallow Deletion",
            "Deletion",
            "Deep Deletion",
        ]
        cnv_cmap = {}
        cnv_names = []
        for name, (key, color) in zip(_cnv_names, _cnv_cmap.items()):
            if key in np.unique(data.cnv.df):
                cnv_cmap[key] = color
                cnv_names.append(name)
        return Palette(cnv_cmap), cnv_names

    def get_tmb_cmap(self, data):
        tmb_cmap = {
            better_effect_legend(effect): self.get(effect, self.grey)
            for effect in reversed(data.tmb.columns)
        } if data.tmb is not None else {}
        return Palette(tmb_cmap)

    def get_mutsig_cmap(self, data):
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
            "clock-like": self.brown,
            "APOBEC": self.red,
            "MMR": self.green,
            "PolE/D/H": self.cyan,
            "Other": self.violet,
            "UV": self.brightorange,
            "Smoking": self.brightblue,
            "Treatment": self.pink,
            "Error": self.grey
        }
        signature_colors = {color: signature_sets[key] for key, color in mutsigset_palette.items()}
        for color, signatures in signature_colors.items():
            for s, c in zip(signatures, self.make_diverging_palette(color=color, n_colors=len(signatures))):
                mutsigset_palette[s] = c
        mutsig_cmap = {
            sig: mutsigset_palette.get(sig, color)
            for sig, color in zip(data.mutsig.columns, sns.color_palette("husl", n_colors=len(data.mutsig.columns)))
        } if data.mutsig is not None else {}
        return Palette(mutsig_cmap)

    def get_meta_cmaps(self, data):
        meta_cmaps = {}

        def add_cmap(col, _palette=None, order=None):
            if col not in data.meta.rows:
                return None

            values = []
            for value in data.meta.df[col].values:
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
            if col not in data.meta.rows:
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
            add_cmap(col, _palette=[self.normalizeRGB(105, 200, 219), self.lightgrey, self.backgroundgray], order=["yes", "no", "unknown"])
        add_cmap("WGD", _palette=self.make_diverging_palette(self.violet, n_colors=5*3)[::5], order=list(range(3)))
        add_cont_cmap("Contamination", plt.cm.get_cmap("BuPu"), 0.0, 0.05)
        add_cont_cmap("Tumor Purity", plt.cm.get_cmap("plasma_r"), 0, 1)
        add_cont_cmap("Ploidy", plt.cm.get_cmap("PiYG"), 1, 6, center=2)
        add_cont_cmap("Subclonal Fraction", plt.cm.get_cmap("RdPu"), 0, 0.5)

        # For all custom columns:
        for col in data.meta.rows:
            if col not in meta_cmaps.keys():
                add_cmap(col)

        return {k: meta_cmaps[k] for k in data.meta.rows if k in meta_cmaps.keys()}

    @staticmethod
    def condense(cmaps):
        continuous_cmaps = {k: p for k, p in cmaps.items() if not isinstance(p, Palette)}
        condenseable_cmaps = {k: p for k, p in cmaps.items() if isinstance(p, Palette)}
        grouped_palettes = defaultdict(list)
        for k, p in condenseable_cmaps.items():
            grouped_palettes[p.hash()].append(k)

        condensed_cmaps = {"\n".join(k_list): Palette.from_hash(h) for h, k_list in grouped_palettes.items()}
        return continuous_cmaps | condensed_cmaps

import numpy as np
import pandas as pd

from comutplotlib.gistic import Gistic
from comutplotlib.maf import MAF
from comutplotlib.sif import SIF


np.random.seed(21)


def make_sif():
    columns = [SIF.sample, SIF.patient, SIF.sample_type, SIF.sex, SIF.material, SIF.histology, SIF.contamination, SIF.tumor_purity, "Paired", SIF.tmb, SIF.tmb_error]
    entries = []
    n_patients = 100
    for i in range(1, n_patients + 1):
        histology = np.random.choice(["BRCA", "LUAD", "MEL", "RCC"], p=[0.4, 0.2, 0.2, 0.2])
        sex = np.random.choice(["Male", "Female", "unknown"], p=[0.4, 0.5, 0.1])
        n_samples = np.random.choice([1, 2, 3, 4], p=[0.8, 0.18, 0.015, 0.005])
        Paired = np.random.choice([True, False], p=[0.7, 0.3])
        for j in range(1, n_samples + 1):
            sample_type = np.random.choice(["BM", "EM"], p=[0.7, 0.3])
            material = np.random.choice(["FFPE", "FF", "cfDNA"], p=[0.5, 0.3, 0.2])
            contamination = np.random.beta(i + 1, n_patients * 10)
            tumor_purity = np.random.beta(3, 3)
            tmb = 19 ** (i / n_patients) + 0.1 * j
            tmb_error = np.abs(np.random.normal(0, 0.1 * tmb))
            entries.append([f"Sample {i}.{j}", f"Patient {i}", sample_type, sex, material, histology, contamination, tumor_purity, Paired, tmb, tmb_error])

    sif = SIF(data=pd.DataFrame(entries, columns=columns))
    sif.to_csv("test.sif.tsv")
    return sif


def make_genes():
    n_genes = 25
    genes = [f"Gene {i}" for i in range(1, n_genes + 1)]
    return genes


def make_maf(sif, genes):
    columns = [MAF.sample, MAF.patient, MAF.chromosome, MAF.start_pos, MAF.end_pos, MAF.gene_name, MAF.type, MAF.ref_allele, MAF.alt_allele, MAF.effect, MAF.ref_count, MAF.alt_count, MAF.protein_change]
    entries = []
    # doesn't matter for comut plot!
    start_pos = 1e6
    end_pos = 1e6 + 1
    ref_count = 1e2
    alt_count = 1e2
    for sample in sif.samples:
        patient = sif.select({SIF.sample: sample}).patients[0]
        for gene in genes:
            id = int(gene.split(" ")[1])
            rate = -1.5 if id in range(1, 13) else -2.1
            for i in range(int(np.random.lognormal(rate, 1))):
                chromosome = gene.split(" ")[1]
                ref_allele, alt_allele = np.random.choice(["A", "C", "G", "T"], size=2, replace=False)
                effect_p = np.arange(10, 1, -1) ** 3
                effect = np.random.choice([MAF.missense, MAF.nonsense, MAF.splice_site, MAF.frame_shift_del, MAF.frame_shift_ins, MAF.utr5, MAF.utr3, MAF.in_frame_del, MAF.in_frame_ins], p=effect_p / effect_p.sum())
                amino_acids = np.random.choice(["A", "C", "D", "E"], size=2, replace=False)
                protein_change = "p." + amino_acids[0] + "42" + amino_acids[1]
                entries.append([sample, patient, chromosome, start_pos, end_pos, gene, "SNV", ref_allele, alt_allele, effect, ref_count, alt_count, protein_change])
            if id == 8 and np.random.rand() < 0.1:
                entries.append([sample, patient, id, start_pos, end_pos, gene, "SNV", "A", "T", MAF.splice_site, ref_count, alt_count, "p.S42P"])
            if id == 10 and np.random.rand() < 0.1:
                entries.append([sample, patient, id, start_pos, end_pos, gene, "SNV", "C", "T", MAF.missense, ref_count, alt_count, "p.M42S"])

    maf = MAF(data=pd.DataFrame(entries, columns=columns))
    maf.to_csv("test.maf.tsv")
    return maf


def make_gistic(sif, genes):
    n_genes = len(genes)
    columns = [Gistic._locus_id, Gistic._cytoband]
    genes = pd.Index(genes, name=Gistic._gene_symbol)
    locus_ids = [f"Locus {i}" for i in range(1, n_genes + 1)]
    cytobands = [f"{i}p1.{i}" for i in range(1, n_genes + 1)]
    thresholded_gistic_scores = [locus_ids, cytobands]
    for patient in sif.patients:
        columns.append(patient)
        neutral_genes = np.random.choice(range(-2, 3), p=[0.001, 0.1, 0.798, 0.1, 0.001], size=8)
        amp_genes = np.random.choice(range(-2, 3), p=[0.01, 0.1, 0.49, 0.28, 0.12], size=6)
        del_genes = np.random.choice(range(-2, 3), p=[0.1, 0.3, 0.49, 0.1, 0.01], size=6)
        other_neutral_genes = np.random.choice(range(-2, 3), p=[0.001, 0.1, 0.798, 0.1, 0.001], size=n_genes - 20)
        thresholded_gistic_scores.append(np.concatenate([neutral_genes, amp_genes, del_genes, other_neutral_genes]))

    gistic = Gistic(data=pd.DataFrame(thresholded_gistic_scores, columns=genes, index=columns).T)
    gistic.data.to_csv("test.all_thresholded.by_genes.txt", sep="\t")
    return gistic


def make_mutsig(sif):
    signatures = ["Ageing", "Smoking", "UV", "Chemotherapy"]
    entries = []
    for patient in sif.patients:
        histology = sif.select({SIF.patient: patient}).data[SIF.histology].unique()[0]
        if histology == "BRCA":
            exposures = np.random.gamma([1, 1, 1, 1], [5, 1, 1, 1])
        elif histology == "LUAD":
            exposures = np.random.gamma([1, 1, 1, 1], [1, 20, 1, 1])
        elif histology == "MEL":
            exposures = np.random.gamma([1, 1, 1, 1], [1, 1, 50, 1])
        elif histology == "RCC":
            exposures = np.random.gamma([1, 1, 1, 1], [1, 1, 1, 2])
        else:
            exposures = np.random.gamma([1, 1, 1, 1])
        entries.append(exposures)

    mutsig = pd.DataFrame(entries, columns=signatures, index=sif.patients)
    mutsig.to_csv("test.mutsig.tsv", sep="\t")
    return mutsig


if __name__ == "__main__":
    sif = make_sif()
    genes = make_genes()
    maf = make_maf(sif, genes)
    gistic = make_gistic(sif, genes)
    mutsig = make_mutsig(sif)

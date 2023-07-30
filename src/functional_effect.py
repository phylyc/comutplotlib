from typing import Iterable
from src.mutation_annotation import MutationAnnotation as MutA


def better_effect_legend(effect):
    return (
        effect
        .replace("_", " ")
        .replace("OutOfFrame", "ooF")
        .replace("InFrame", "iF")
    )


def sort_functional_effects(effects: Iterable[str], ascending: bool = True) -> list[str]:
    """ Sort a list of functional effects based on a gain of function - loss of function scale.
    :param effects:
    :param ascending: Sort loss of function first, then towards gain of function.
    :return: Sorted list of functional effects.
    """

    function_values = {
        MutA.gain_of_function: 1000,
        MutA.missense: 100,
        MutA.in_frame_ins: 52,
        MutA.in_frame_del: 51,
        MutA.de_novo_start_in_frame: 20,

        MutA.utr5: 11,
        MutA.utr3: 10,
        MutA.silent: 5,
        MutA.synonymous: 0,
        MutA.intron: 0,
        MutA.flank5: -3,
        MutA.flank3: -4,
        MutA.igr: -5,

        MutA.rna: -30,
        MutA.linc_rna: -40,

        MutA.frame_shift_ins: -51,
        MutA.frame_shift_del: -52,
        MutA.de_novo_start_out_of_frame: -59,
        MutA.start_codon_snp: -60,
        MutA.start_codon_ins: -61,
        MutA.start_codon_del: -62,
        MutA.translation_start_site: -70,
        MutA.splice_site: -80,
        MutA.nonstop: -90,
        MutA.read_through: -91,
        MutA.stop_codon_del: -95,
        MutA.nonsense: -100,
        MutA.stop_codon_ins: -105,
        MutA.nonsynonymous: -666,
        MutA.structural: -999,
        MutA.loss_of_function: -1000,

        MutA.snv: 1,
        MutA.dnv: 2,
        MutA.tnv: 3,
        MutA.mnv: 4,
        MutA.insertion: 5,
        MutA.deletion: 6
    }

    return list(sorted(effects, key=lambda e: function_values.get(e, 0), reverse=ascending))

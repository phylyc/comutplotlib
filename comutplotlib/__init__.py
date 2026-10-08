from comutplotlib.comut_argparse import parse_args
from comutplotlib.comut import Comut
from comutplotlib.comut_data import ComutData
from comutplotlib.comut_layout import ComutLayout
from comutplotlib.comut_panels import ComutPanels, DEFAULT_PANELS
from comutplotlib.comut_plotter import ComutPlotter
from comutplotlib.functional_effect import sort_functional_effects

from comutplotlib.gistic import Gistic, join_gistics
from comutplotlib.mark import Mark, join_marks
from comutplotlib.seg import SEG, join_segs
from comutplotlib.cnv import CNV
from comutplotlib.epi import EPI

from comutplotlib.maf import MAF, join_mafs
from comutplotlib.maf_encoding import MAFEncoding
from comutplotlib.mutation_annotation import MutationAnnotation
from comutplotlib.mutational_signature_set import MutationalSignatureSet
from comutplotlib.snv import SNV

from comutplotlib.sif import SIF, join_sifs
from comutplotlib.sample_annotation import SampleAnnotation
from comutplotlib.meta import Meta

from comutplotlib.annotation_table import AnnotationTable
from comutplotlib.layout import Layout
from comutplotlib.panel import Panel
from comutplotlib.palette import Palette
from comutplotlib.plotter import Plotter

from comutplotlib.mathutils import decompose_rectangle_into_polygons
from comutplotlib.pandas_util import *

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("comutplotlib")
except PackageNotFoundError:
    # Package metadata is unavailable when running from a source checkout
    # that has not been installed (e.g. `pip install -e .`).
    __version__ = "0.0.0+unknown"
__all__ = ["__version__"]

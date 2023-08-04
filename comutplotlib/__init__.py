from comutplotlib.comut_argparse import parse_args
from comutplotlib.comut import Comut
from comutplotlib.comut_data import ComutData
from comutplotlib.comut_layout import ComutLayout
from comutplotlib.comut_plotter import ComutPlotter
from comutplotlib.functional_effect import sort_functional_effects

from comutplotlib.gistic import Gistic, join_gistics
from comutplotlib.seg import SEG, join_segs
from comutplotlib.cnv import CNV

from comutplotlib.maf import MAF, join_mafs
from comutplotlib.maf_encoding import MAFEncoding
from comutplotlib.mutation_annotation import MutationAnnotation
from comutplotlib.snv import SNV

from comutplotlib.sif import SIF, join_sifs
from comutplotlib.sample_annotation import SampleAnnotation
from comutplotlib.meta import Meta

from comutplotlib.annotation_table import AnnotationTable
from comutplotlib.layout import Layout
from comutplotlib.panel import Panel
from comutplotlib.palette import Palette
from comutplotlib.plotter import Plotter

from comutplotlib.math import decompose_rectangle_into_polygons
from comutplotlib.pandas_util import *

from pkg_resources import get_distribution

__version__ = get_distribution("comutplotlib").version
__all__ = ["__version__"]

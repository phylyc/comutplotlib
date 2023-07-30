from src.comut_argparse import parse_args
from src.comut import Comut


args = parse_args()
comut = Comut(**vars(args))
comut.make_comut()

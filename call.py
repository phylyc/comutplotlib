from comutplotlib.comut_argparse import parse_args
from comutplotlib.comut import Comut


def main():
    args = parse_args()
    comut = Comut(**vars(args))
    comut.make_comut()


if __name__ == "__main__":
    main()


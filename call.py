from src.comut_argparse import parse_args
from src.plot_comutation_table import plot_comutation_table


args = parse_args()

print("args:")
for i, arg in enumerate(args._get_args()):
    print(f"{i}\t{arg}")
print("kwargs:")
for i, (kw, arg) in enumerate(vars(args).items()):
    print(f"{i}\t{kw}: {arg}")


plot_comutation_table(**vars(args))

import re


class AminoAcidChange(object):
    @staticmethod
    def from_string(string: str) -> "AminoAcidChange":
        """string has the pattern 'p.NxxxM', where N is the code of the
        base amino acid, M is the code of the changed amino acid, and
        xxx is the index position of the change in the chain.
        """
        array = re.split("(\d+)", string.split(".")[1])
        base = array[0]
        pos = array[1]
        change = array[2]
        return AminoAcidChange(base=base, change=change, pos=pos)

    def __init__(self, base: str, change: str, pos=None) -> None:
        self.base = AminoAcid(code=base)
        self.change = AminoAcid(code=change)
        self.pos = pos


class AminoAcid(object):
    def __init__(self, code: str) -> None:
        self.code = code

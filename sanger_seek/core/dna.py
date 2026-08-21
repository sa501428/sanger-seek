"""Nucleotide / amino-acid utilities."""

from __future__ import annotations

_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "D": "H", "H": "D", "V": "B", "N": "N", "-": "-",
}

IUPAC_SETS: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

SETS_TO_IUPAC: dict[frozenset[str], str] = {
    frozenset(v): k for k, v in IUPAC_SETS.items()
}

VALID_BASES = set(IUPAC_SETS) | {"-"}

# Extra equalities so ambiguity codes in reads match plain bases in the
# reference during alignment (edlib additionalEqualities format).
AMBIG_EQUALITIES: list[tuple[str, str]] = [
    (code, base)
    for code, bases in IUPAC_SETS.items()
    if len(bases) > 1
    for base in bases
]


def complement(base: str) -> str:
    return _COMPLEMENT.get(base.upper(), "N")


def revcomp(seq: str) -> str:
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def expand(base: str) -> str:
    """IUPAC code -> string of concrete bases it can represent."""
    return IUPAC_SETS.get(base.upper(), "")


def bases_match(a: str, b: str) -> bool:
    """True if the IUPAC sets of two base codes intersect."""
    sa, sb = expand(a), expand(b)
    return bool(set(sa) & set(sb))


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA_THREE = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "X": "Xaa",
}


def translate_codon(codon: str) -> str:
    return CODON_TABLE.get(codon.upper(), "X")


def translate(seq: str) -> str:
    return "".join(
        translate_codon(seq[i : i + 3]) for i in range(0, len(seq) - 2, 3)
    )


def aa3(aa1: str) -> str:
    return AA_THREE.get(aa1.upper(), "Xaa")

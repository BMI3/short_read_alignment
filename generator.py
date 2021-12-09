import random


nucleotides = ("A", "T", "C", "G")


def generate_DNA(length):
    return "".join(map(lambda x: nucleotides[x], random.choices(range(0, 4), k=length)))

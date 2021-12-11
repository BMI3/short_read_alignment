import random

nucleotides = ("A", "T", "C", "G")


def generate_DNA(length):
    return "".join(map(lambda x: nucleotides[x], random.choices(range(0, 4), k=length)))


def extract(seq, lm, n):
    """randomly extract n motifs with length of lm from the reference"""
    index0 = random.sample(range(len(seq)), n)
    motifs = [seq[i : min(i + lm, len(seq))] for i in index0]
    return motifs, index0


def mutate(s, n):
    """
    randomly mutate n position in s and return mutated s without changing s.
    if n >= len(s), all the position would be mutate.
    mutation uses only characters appeared in s.
    """
    characters = set(s)
    n = min(n, len(s))
    position = random.sample(range(len(s)), n)
    ms = s[:]
    for p in position:
        ms[p] = random.choice(
            characters
            - {
                s[p],
            }
        )
    return ms

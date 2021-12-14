from MyBWT import MyBWT
import pickle
import sys
import argparse

fa_set = {"fa", "fasta"}
fq_set = {"fq", "fastq"}


def readGenome(filename):
    """
    read fa file
    """
    genome = ""
    with open(filename, "r") as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == ">":
                genome += line.rstrip()
    my_gonome = MyBWT(genome)
    with open(filename + ".bwt", "wb") as f:
        pickle.dump(my_gonome, f, 0)
    return my_gonome


def reLoadRefObj(filename):
    with open(filename, mode="rb") as f:
        my_gonome = pickle.load(f, encoding="UTF-8")

    print("debug here to check the object")

    print(sys.getsizeof(my_gonome.lc))
    print(sys.getsizeof("".join(my_gonome.lc)))


parser = argparse.ArgumentParser()
parser.add_argument(
    "-b",
    "--build",
    action="store_true",
    help="build a MyBWT object of the input reference genome",
)
parser.add_argument(
    "-q",
    "--query",
    action="store_true",
    help="align short reads to the genome",
)
parser.add_argument(
    "-r",
    "--reference",
    required=True,
    help="input reference genome (.fa or .fa.bwt)",
)
parser.add_argument(
    "-s",
    "--shortRead",
    help="input sequncing reads (.fq)",
)

if __name__ == "__main__":
    args = parser.parse_args()
    fns = args.reference.split(".")
    if args.build:
        if fns[-1] in fa_set:
            print("start")
            readGenome(args.reference)
            print("finished")
        else:
            print("input file should be fasta file")
    elif args.query and args.shortRead:
        if fns[-1] == "bwt" and fns[-2] in fa_set:
            genome = reLoadRefObj(args.reference)
        elif fns[-1] in fa_set:
            genome = readGenome(args.reference)
        else:
            print("unrecognized genome file")
        genome.seeding

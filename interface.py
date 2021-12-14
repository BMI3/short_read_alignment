from MyBWT import MyBWT
import pickle
import sys
import argparse


def readGenome(filename):
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
    help="input reference genome (.fa or .fa.bwt)",
)

if __name__ == "__main__":
    args = parser.parse_args()
    if args.build:
        print("start")
        readGenome(args.reference)
        print("finished")
    elif args.query and args.shortRead:
        if (
            args.reference.split(".")[-1] == "bwt"
            and args.reference.split(".")[-2] == "fa"
        ):
            genome = reLoadRefObj(args.reference)
        elif args.reference.split(".")[-1] == "fa":
            genome = readGenome(args.reference)
        else:
            print("unrecognized genome file")
        genome.seeding

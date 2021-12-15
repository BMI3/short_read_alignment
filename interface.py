from MyBWT import MyBWT
import pickle
import time
import argparse

# from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

fa_set = {"fa", "fasta"}
fq_set = {"fq", "fastq"}


def readGenome(filename):
    """
    read fa file
    """
    time0 = time.time()
    print("start")
    # genome = ""
    # with open(filename, "r") as f:
    #     for line in f:
    #         # ignore header line with genome information
    #         if not line[0] == ">":
    #             genome += line.rstrip()
    # my_gonome = MyBWT(genome)
    my_gonome = MyBWT(str(next(SeqIO.parse(filename, "fasta")).seq))
    with open(filename + ".bwt", "wb") as f:
        pickle.dump(my_gonome, f, 0)
    print("finish buiding the genome index, cost", time.time() - time0, "s")
    return my_gonome


def reLoadRefObj(filename):
    with open(filename, mode="rb") as f:
        my_gonome = pickle.load(f, encoding="UTF-8")

    # print("debug here to check the object")
    # print(sys.getsizeof(my_gonome.lc))
    # print(sys.getsizeof("".join(my_gonome.lc)))
    return my_gonome


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
    "-p",
    "--pairend",
    action="store_true",
    help="if the sequencing data is pair-end",
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

            readGenome(args.reference)

        else:
            print("input should be fasta file")
    elif args.query and args.shortRead:
        if fns[-1] == "bwt":
            genome = reLoadRefObj(args.reference)
        elif fns[-1] in fa_set:
            genome = readGenome(args.reference)
        else:
            print("unrecognized genome file")

        # with open(args.shortRead) as in_handle:
        #     for title, seq, qual in FastqGeneralIterator(in_handle):
        #         genome.seeding(seq)

        fns = args.shortRead.split(".")[-1]
        if fns in fa_set:
            suffix = "fasta"
        if fns in fq_set:
            suffix = "fastq"
        records = list(SeqIO.parse(args.shortRead, suffix))
        print(len(records))
        no = 0
        for i, record in enumerate(records):
            print(record.name)
            result = genome.seeding(record.seq)
            if args.pairend:
                result += genome.seeding(record.reverse_complement().seq)
            if result == []:
                no += 1
            print("result:", result)
            print("----------------------")
        print(no, "/", len(records))

from MyBWT import MyBWT
import pickle


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


def reLoadRefObj(filename):
    with open(filename, mode="rb") as f:
        my_gonome = pickle.load(f, encoding="UTF-8")

    print("debug here to check the object")


# readGenome("./data/phix.fa")
reLoadRefObj("./data/phix.fa.bwt")

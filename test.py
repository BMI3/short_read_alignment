import time
import os
from numpy import cos
import psutil
import matplotlib.pyplot as plt
import ag2
import ag_tra
import generator
from bwtclass import BWTtrans


def ag2_vs_agt(iter, show=False, mode="t"):
    cost = []
    for i in iter:
        print(i)
        s = generator.generate_DNA(i)
        if mode == "t":
            time0 = time.time()
            ag2.BWTstore(s)
            time1 = time.time()
            ag_tra.fl(s)
            time2 = time.time()
            c = (time1 - time0, time2 - time1)
        else:
            pid = os.getpid()
            p = psutil.Process(pid)
            memory0 = p.memory_full_info().uss / 1024
            ag2.BWTstore(s)
            memory1 = p.memory_full_info().uss / 1024
            ag_tra.fl(s)
            memory2 = p.memory_full_info().uss / 1024
            c = (memory1 - memory0, memory2 - memory1)
        cost.append(c)
        print(c)
    nn = len(cost)
    plt.xscale("log")
    plt.plot(
        iter,
        [cost[i][0] for i in range(nn)],
        iter,
        [cost[i][1] for i in range(nn)],
    )
    if show:
        plt.show()
    else:
        plt.savefig(str(iter) + ".png")


def test_lc():
    """test our function to get the last column of the sorted N*N matrix with simply generate the whole matrix, sort, and extract"""
    logs = [10 ** i for i in range(2, 5)]
    ag2_vs_agt(logs, mode="m")

    ag2_vs_agt(range(10 ** 4, 7 * 10 ** 4, 10 ** 4))
    ag2_vs_agt(range(10 ** 4, 2 * 10 ** 4, 10 ** 3))
    ag2_vs_agt(range(16 * 10 ** 3, 17 * 10 ** 3, 10 ** 2))
    ag2_vs_agt(range(163 * 10 ** 2, 164 * 10 ** 2, 10))


def find_vs_query(obj, motif):
    time0 = time.time()
    obj[0].querySA(motif[0][0])
    time1 = time.time()
    obj[1].find(motif)
    time2 = time.time()
    return (time1 - time0, time2 - time1)


def find_vs_query_memory(s, motif):
    pid = os.getpid()
    p = psutil.Process(pid)
    memory0 = p.memory_full_info().uss / 1024
    BWTtrans(s, motif).querySA(motif[0][0])
    memory1 = p.memory_full_info().uss / 1024
    ag2.BWTstore(s).find(motif)
    memory2 = p.memory_full_info().uss / 1024
    return (memory1 - memory0, memory2 - memory1)


def test_query(mm, mode="t"):
    for i in range(2, 6):
        print(i)
        s = generator.generate_DNA(10 ** i)
        m = generator.extract(s, 10 ** (i - 1), mm)
        obj = BWTtrans(s, ""), ag2.BWTstore(s)

        cost = []
        for j in range(mm):
            if mode == "t":
                c = find_vs_query(obj, m[0][j])
            else:
                c = find_vs_query_memory(s, m[0][j])
            cost.append(c)
            print(c)


if __name__ == "__main__":
    # test_lc()
    # test_query(5, mode="m")
    test_query(5, mode="t")

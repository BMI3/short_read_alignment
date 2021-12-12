import time
import os
import psutil
import matplotlib.pyplot as plt
import ag2
import ag_tra
import generator


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


if __name__ == "__main__":
    test_function()

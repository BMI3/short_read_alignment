import time
import matplotlib.pyplot as plt
import sys
import logging
import ag2
import ag_tra
import generator

# logging.basicConfig


"""decerator"""


def evalue(function):
    def timer(*args, **kwargs):
        time0 = time.time()
        result = function(*args, **kwargs)
        time_cost = time.time() - time0
        print("time cost:", time_cost)
        return result

    def memory(*args, **kwargs):
        import os
        import psutil

        pid = os.getpid()
        p = psutil.Process(pid)
        info_start = p.memory_full_info().uss / 1024
        result = function(*args, **kwargs)
        info_end = p.memory_full_info().uss / 1024
        print(str(info_end - info_start) + "KB")
        return result

    # if mode == "t":
    #     return timer
    # if mode == "m":
    #     return memory
    return timer


@evalue
def ag2_fl(s):
    print("ag2")
    ag2.BWTstore(s)


@evalue
def agt(s):
    print("agt")
    ag_tra.fl(s)


def ag2_vs_agt(iter, show=False):
    time_cost = []
    for i in iter:
        print(i)
        s = generator.generate_DNA(i)
        time0 = time.time()
        ag2.BWTstore(s)
        time1 = time.time()
        ag_tra.fl(s)
        time2 = time.time()
        tc = time1 - time0, time2 - time1
        time_cost.append(tc)
        print(tc)
    nn = len(time_cost)
    plt.xscale("log")
    plt.plot(
        iter,
        [time_cost[i][0] for i in range(nn)],
        iter,
        [time_cost[i][1] for i in range(nn)],
    )
    if show:
        plt.show()
    else:
        plt.savefig(str(iter) + ".png")


def test_function():
    logs = [10 ** i for i in range(2, 8)]
    ag2_vs_agt(logs)

    ag2_vs_agt(range(10 ** 4, 8 * 10 ** 4, 10 ** 4))
    ag2_vs_agt(range(10 ** 4, 2 * 10 ** 4, 10 ** 3))
    ag2_vs_agt(range(16 * 10 ** 3, 17 * 10 ** 3, 10 ** 2))
    ag2_vs_agt(range(163 * 10 ** 2, 164 * 10 ** 2, 10))


if __name__ == "__main__":
    test_function()
    # print(sys.platform == "win32")

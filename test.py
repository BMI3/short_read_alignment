import time
import matplotlib.pyplot as plt
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


def test_function():
    time_cost = []
    try:
        for i in range(2, 5):
            print(i)
            s = generator.generate_DNA(10 ** i)
            time0 = time.time()
            ag2.BWTstore(s)
            time1 = time.time()
            ag_tra.fl(s)
            time2 = time.time()
            time_cost.append((time1 - time0, time2 - time1))
            print(time1 - time0, time2 - time1)
    except MemoryError:
        print("memoryError at i=", i)
    # fig, ax = plt.subplots(1, 2)
    # ax[0].plot(range(2, 5), [time_cost[i][0] for i in range(len(time_cost))])
    # ax[1].plot(range(2, 5), [time_cost[i][1] for i in range(len(time_cost))])
    nn = len(time_cost)
    plt.plot(
        range(2, 2 + nn),
        [time_cost[i][0] for i in range(nn)],
        range(2, 2 + nn),
        [time_cost[i][1] for i in range(nn)],
    )
    plt.show()


if __name__ == "__main__":
    test_function()

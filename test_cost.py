import time
import os
from numpy.lib.function_base import append
import psutil
import numpy as np
import matplotlib.pyplot as plt
from MyBWT import MyBWT
import test_generator


def evalue(mode="t"):
    def timer(f):
        def wrapper(*args, **kwargs):
            t0 = time.time()
            f(*args, **kwargs)
            t1 = time.time()
            return t1 - t0

        return wrapper

    def memory(f):
        def wrapper(*args, **kwargs):
            pid = os.getpid()
            p = psutil.Process(pid)
            m0 = p.memory_full_info().uss / 1024
            f(*args, **kwargs)
            m1 = p.memory_full_info().uss / 1024
            return m1 - m0

        return wrapper

    if mode == "t":
        return timer
    else:
        return memory


def cost_generate_LC(iter, mode="t", logs=True, fs=6):
    @evalue(mode)
    def generate_LC_rotation(s):
        s += "$"
        result = list(map(lambda i: s[-1 - i :] + s[: -1 - i], range(len(s))))
        return "".join([x[-1] for x in sorted(s)])

    @evalue(mode)
    def generate_LC_SA(s):
        s += "$"
        sa = []
        for i in range(len(s)):
            sa.append((s[i:], i))
        lc = ""
        for i in sa:
            if i[1] == 0:
                lc += "$"
            else:
                lc += s[i[1] - 1]
        return lc

    @evalue(mode)
    def generate_LC(s):
        MyBWT.generateLC(MyBWT("", ""), s)

    rotation = []
    SA = []
    my = []
    for i in iter:
        print(i)
        s = test_generator.generate_DNA(i)
        if fs % 2:
            rotation.append(generate_LC_rotation(s))
        if fs % 4 // 2:
            SA.append(generate_LC_SA(s))
        my.append(generate_LC(s))

    return my, SA, rotation


def graph_generate_LC(iter, data, show=False, mode="t", logs=True, fs=6):
    ax = plt.subplots()

    if logs:
        x = np.arange(len(iter))
        width = 0.3
        ax.set_yscale("log")
        if fs % 2:
            ax1 = ax.bar(x - width, data[2], width, label="rotation", color="#F7D6AA")
        if fs % 4 // 2:
            ax2 = ax.bar(x + 1 - 1, data[1], width, label="SA", color="#4D454D")
        ax3 = ax.bar(x + width, data[0], width, label="my", color="#C95B45")
        ax.set_xticks(x)
        ax.set_xticklabels(iter)
    else:
        if fs % 2:
            ax1 = ax.plot(iter, data[2], label="rotation", color="#F7D6AA")
        if fs % 4 // 2:
            ax2 = ax.plot(iter, data[1], label="SA", color="#4D454D")
        ax3 = ax.plot(iter, data[0], label="my", color="#C95B45")

    ax.set_xlabel("length of reference sequence")
    if mode == "t":
        ax.set(
            ylabel="time cost/s",
            title="time cost of methods to construct BWT",
        )
    else:
        ax.set(
            ylabel="memory cost/KB",
            title="memory cost of methods to construct BWT",
        )

    ax.legend()

    # ax.bar_label(ax1, padding=3)
    # ax.bar_label(ax2, padding=3)
    # ax.bar_label(ax3, padding=3)

    if show:
        plt.show()
    else:
        fn = "b" + mode + str(iter[0]) + ".png"
        print("save as", fn)
        plt.savefig(fn)


def test_generate_LC():
    """
    compare our algorithm with
        1. generate with the whole rotation matrix
        2. generate with a suffix tree
    """

    # logs = [10 ** i for i in range(2, 6)]
    # time_logs = cost_generate_LC(logs, mode="t", fs=4)
    # graph_generate_LC(logs, time_logs, mode="t", fs=4)
    # memory_logs = cost_generate_LC(logs, mode="m", fs=4)
    # graph_generate_LC(logs, memory_logs, mode="m", fs=4)

    def repeat(gradient, n, mode="t", fs=6):
        choose = fs // 4 + fs % 4 // 2 + fs % 2
        for i in range(n):
            if i == 0:
                result = np.array(cost_generate_LC(gradient, mode, logs=False)[:choose])
            else:
                result += np.array(
                    cost_generate_LC(gradient, mode, logs=False)[:choose]
                )
        return result / n

    gradient = range(6 * 10 ** 4, 7 * 10 ** 4, 5 * 10 ** 2)
    time_gradient = repeat(gradient, 3)
    graph_generate_LC(gradient, time_gradient, mode="t", logs=False)
    memory_gradient = repeat(gradient, 3, mode="m")
    graph_generate_LC(gradient, memory_gradient, mode="m", logs=False)


def cost_query(iter, mode="t"):
    """
    compare query with tally with
            query with rank
    """

    @evalue(mode)
    def query_with_rank(bwt, sr_list):
        for sr in sr_list:
            result = []
            possible = range(bwt.char_range[sr[-1]][0], bwt.char_range[sr[-1]][1])
            for i in possible:
                next = i
                go = 2
                while len(sr) > go and bwt.lc[next] == sr[-go]:
                    chara = bwt.lc[next]
                    rank = bwt.tally[chara][next]
                    next = bwt.char_range[chara][0] + rank - 1
                    go += 1
                if go == len(sr):
                    result.append(bwt.po[i] - len(sr) + 1)

    @evalue(mode)
    def query_with_tally(bwt, sr_list):
        for sr in sr_list:
            bwt.query(sr)

    ref = test_generator.generate_DNA(10000)
    bwt = MyBWT("", ref)
    rank = []
    tally = []
    for i in iter:
        print(i)
        sr_list = test_generator.extract(ref, 50, i)[0]
        rank.append(query_with_rank(bwt, sr_list))
        tally.append(query_with_tally(bwt, sr_list))
    return rank, tally


def graph_query(iter, data, show=False, mode="t"):
    fig, ax = plt.subplots()
    ax.plot(iter, data[0], label="rank", color="#4D454D")
    ax.plot(iter, data[1], label="tally", color="#C95B45")
    ax.set_xlabel("length of reference sequence")
    if mode == "t":
        ax.set(
            ylabel="time cost/s",
            title="time cost of methods to query",
        )
    else:
        ax.set(
            ylabel="memory cost/KB",
            title="memory cost of methods to query",
        )
    ax.legend()

    if show:
        plt.show()
    else:
        fn = "q" + mode + str(iter[0]) + ".png"
        print("save as", fn)
        plt.savefig(fn)


def test_query():
    iter = range(5 * 10 ** 2, 5 * 10 ** 3, 5 * 10 ** 2)
    time_cost = cost_query(iter)
    graph_query(iter, time_cost)


if __name__ == "__main__":
    # test_generate_LC()
    test_query()

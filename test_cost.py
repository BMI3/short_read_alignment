import time
import os
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


def cost_generate_LC(iter, mode="t", logs=True):
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
        MyBWT.generateLC(MyBWT(""), s)

    rotation = []
    SA = []
    my = []
    for i in iter:
        print(i)
        s = test_generator.generate_DNA(i)
        # if logs and mode == "t":
        #     rotation.append(generate_LC_rotation(s))

        SA.append(generate_LC_SA(s))
        my.append(generate_LC(s))

    return SA, my, rotation


def graph_generate_LC(iter, data, show=False, mode="t", logs=True):
    fig, ax = plt.subplots()

    if logs:
        x = np.arange(len(iter))
        width = 0.3
        ax.set_yscale("log")
        if mode == "t":
            ax1 = ax.bar(x - width, data[2], width, label="rotation", color="#F7D6AA")
        ax2 = ax.bar(x + 1 - 1, data[0], width, label="SA", color="#4D454D")
        ax3 = ax.bar(x + width, data[1], width, label="my", color="#C95B45")
        ax.set_xticks(x)
        ax.set_xticklabels(iter)
    else:
        ax2 = ax.plot(iter, data[0], label="SA", color="#4D454D")
        ax3 = ax.plot(iter, data[1], label="my", color="#C95B45")

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
        fn = mode + str(iter[0]) + ".png"
        print("save as", fn)
        plt.savefig(fn)


def test_generate_LC():
    """
    compare our algorithm with
        1. generate with the whole rotation matrix
        2. generate with a suffix tree
    """

    logs = [10 ** i for i in range(2, 7)]
    time_logs = cost_generate_LC(logs, mode="t")
    graph_generate_LC(logs, time_logs, mode="t")
    memory_logs = cost_generate_LC(logs, mode="m")[:2]
    graph_generate_LC(logs, memory_logs, mode="m")

    def repeat(gradient, n):
        for i in range(n):
            if i == 0:
                result = np.array(cost_generate_LC(gradient, mode="t", logs=False)[:2])
            else:
                result += np.array(cost_generate_LC(gradient, mode="t", logs=False)[:2])
        return result / n

    gradient = range(6 * 10 ** 4, 7 * 10 ** 4, 5 * 10 ** 2)
    time_gradient = repeat(gradient, 5)
    graph_generate_LC(gradient, time_gradient, mode="t", logs=False)


# def find_vs_query(obj, motif):
#     time0 = time.time()
#     obj[0].querySA(motif[0][0])
#     time1 = time.time()
#     obj[1].find(motif)
#     time2 = time.time()
#     return (time1 - time0, time2 - time1)


# def find_vs_query_memory(s, motif):
#     pid = os.getpid()
#     p = psutil.Process(pid)
#     memory0 = p.memory_full_info().uss / 1024
#     BWTtrans(s, motif).querySA(motif[0][0])
#     memory1 = p.memory_full_info().uss / 1024
#     ag2.BWTstore(s).find(motif)
#     memory2 = p.memory_full_info().uss / 1024
#     return (memory1 - memory0, memory2 - memory1)


# def test_query(mm, mode="t"):
#     for i in range(2, 6):
#         print(i)
#         s = generator.generate_DNA(10 ** i)
#         m = generator.extract(s, 10 ** (i - 1), mm)
#         obj = BWTtrans(s, ""), ag2.BWTstore(s)

#         cost = []
#         for j in range(mm):
#             if mode == "t":
#                 c = find_vs_query(obj, m[0][j])
#             else:
#                 c = find_vs_query_memory(s, m[0][j])
#             cost.append(c)
#             print(c)


if __name__ == "__main__":
    test_generate_LC()
    # test_query(5, mode="m")
    # test_query(5, mode="t")

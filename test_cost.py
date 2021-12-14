import time
import os
import psutil
import numpy as np
import matplotlib.pyplot as plt
from MyBWT import MyBWT
import generator


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


def cost_generate_LC(iter, show=False, mode="t"):
    @evalue(mode)
    def generate_LC_matrix(s):
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

    matrix = []
    SA = []
    my = []
    for i in iter:
        print(i)
        s = generator.generate_DNA(i)
        matrix.append(generate_LC_matrix(s))
        SA.append(generate_LC_SA(s))
        my.append(generate_LC(s))

    print(matrix, SA, my)

    x = np.arange(len(iter))
    width = 0.3

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax1 = ax.bar(x - width, matrix, width, label="matrix", color="#F7D6AA")
    ax2 = ax.bar(x, SA, width, label="SA", color="#4D454D")
    ax3 = ax.bar(x + width, my, width, label="my", color="#C95B45")

    ax.set_ylabel("time cost/s")
    ax.set_title("time cost of 3 method to generate the BWT last column")
    ax.set_xticks(x, iter)
    ax.legend()

    ax.bar_label(ax1, padding=3)
    ax.bar_label(ax2, padding=3)
    ax.bar_label(ax3, padding=3)

    fig.tight_layout()
    if show:
        plt.show()
    else:
        plt.savefig(str(iter) + ".png")


def test_generate_LC():
    """
    compare our algorithm with
        1. generate with the whole rotation matrix
        2. generate with a suffix tree
    """

    logs = [10 ** i for i in range(2, 5)]
    cost_generate_LC(logs, mode="m", show=True)

    cost_generate_LC(range(10 ** 4, 7 * 10 ** 4, 10 ** 4), mode="t", show=True)

    # ag2_vs_agt(range(10 ** 4, 7 * 10 ** 4, 10 ** 4))
    # ag2_vs_agt(range(10 ** 4, 2 * 10 ** 4, 10 ** 3))
    # ag2_vs_agt(range(16 * 10 ** 3, 17 * 10 ** 3, 10 ** 2))
    # ag2_vs_agt(range(163 * 10 ** 2, 164 * 10 ** 2, 10))


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

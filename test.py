import ag2
import ag_tra
import generator

"""decerator"""


def evalue(function, mode):
    def timer(*args, **kwargs):
        import time

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

    if mode == "t":
        return timer
    if mode == "m":
        return memory
    return function


# @evalue(mode="t")
def ag2_fl(s):
    print("ag2")
    ag2.BWTstore(s)


# @evalue
def agt(s):
    print("agt")
    ag_tra.fl(s)


# def test_function():
#     s = generator.generate_DNA(1000)
#     ag2_fl(s)
#     agt(s)

#         wrong = 0
#     ee = 0
#     for i in range(1000):
#         a = generator.generate_DNA(1000)
#         # print(a)
#         # print(fl(a))
#         # print(ag_tra.fl(a))

#         try:
#             if "".join(BWTstore(a).lc) != ag_tra.fl(a):
#                 print(a, "".join(BWTstore(a).lc), ag_tra.fl(a))
#                 wrong += 1
#         except:
#             print(a, "   error")
#             ee += 1

#     print(wrong, ee)


if __name__ == "__main__":
    # test_function()
    # ag2_fl(generator.generate_DNA(1000))

    s = generator.generate_DNA(1000)
    ex = generator.extract(s, 100, 2)
    print(ex)
    print(ag2.BWTstore(s).find(ex[0][0]))
    print(ag2.BWTstore(s).find(ex[0][1]))
    print(ag2.BWTstore(s).seeding(ex[0][0], 0))
    print(ag2.BWTstore(s).seeding(ex[0][1], 0))

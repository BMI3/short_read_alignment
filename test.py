try:
    import ag2
    import ag_tra
    import generator
except ModuleNotFoundError:
    from pip._internal import main

    main.main(["install", "-r", "requirements.txt"])


"""decerator"""


def timer(function):
    import time

    def wrapper_c(*args, **kwargs):
        time0 = time.time()
        result = function(*args, **kwargs)
        time_cost = time.time() - time0
        print("time cost:", time_cost)
        return result

    return wrapper_c


def count_info(func):
    import os
    import psutil

    def float_info(*args, **kwargs):
        pid = os.getpid()
        p = psutil.Process(pid)
        info_start = p.memory_full_info().uss / 1024
        func(*args, **kwargs)
        info_end = p.memory_full_info().uss / 1024
        print("程序占用了内存" + str(info_end - info_start) + "KB")

    return float_info


# @timer
@count_info
def ag2_fl(s):
    print("ag2")
    ag2.BWTstore(s)


# @timer
@count_info
def agt(s):
    print("agt")
    ag_tra.fl(s)


@try_module
def test_function():
    s = generator.generate_DNA(100000)
    ag2_fl(s)
    agt(s)


if __name__ == "__main__":
    test_function()

# if __name__ == "__main__":
#     wrong = 0
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

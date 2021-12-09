import ag2
import ag_tra
import generator


def timer(function):
    import time

    def wrapper_c(*args, **kwargs):
        time0 = time.time()
        result = function(*args, **kwargs)
        time_cost = time.time() - time0
        print("time cost:", time_cost)
        return result

    return wrapper_c


@timer
def ag2_fl(s):
    print("ag2")
    ag2.BWTstore(s)


@timer
def agt(s):
    print("agt")
    ag_tra.fl(s)


if __name__ == "__main__":
    s = generator.generate_DNA(100000)
    ag2_fl(s)
    agt(s)

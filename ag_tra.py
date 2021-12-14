import numpy as np
import generator


def rong(s):
    s += "$"
    result = list(map(lambda i: s[-1 - i :] + s[: -1 - i], range(len(s))))
    print("rong", result)
    return result


def sorted_rong(s):
    ll = list(rong(s))
    ll.sort()
    return ll


def fl(s):
    return "".join([x[-1] for x in sorted_rong(s)])


def arg(s):
    plus = [x + str(i) for i, x in enumerate(rong(s))]
    plus.sort()
    return [x[-1] for x in plus]


if __name__ == "__main__":
    c = "ACGATATG"
    r = rong(c)
    for i in range(len(r)):
        print(i, r[i])

import numpy as np


class BWTstore:
    def __init__(self, seq) -> None:
        lc = list("$" + seq)
        fc = list(seq + "$")
        self.po = np.argsort(fc)
        self.lc = np.array(lc)[self.po]

        # print(sorted(fc))
        # print(list(self.lc))

        long = len(self.lc)
        ns = [0] * long
        counter = {}
        first = {}
        for i in range(long - 1, -1, -1):
            if self.lc[i] not in counter.keys():
                counter[self.lc[i]] = 1
            else:
                ns[i] = counter[self.lc[i]]
                counter[self.lc[i]] += 1
        start = 0
        self.characters = sorted(counter.keys())
        for i in self.characters:
            first[i] = start
            start += counter[i]

        # print(list(map(str, ns)))
        self.ns = ns
        self.first = first

    def find(self, motif):
        result = []
        possible = range(self.first[motif[-1]], self.next_first(motif[-1]))
        for i in possible:
            next = i
            go = 2
            while len(motif) > go and self.lc[next] == motif[-go]:
                ns = self.ns[next]
                chara = self.lc[next]
                next = self.next_first(chara) - ns - 1
                go += 1
            if go == len(motif):
                result.append(self.po[i] - len(motif) + 1)
        return result

    def next_first(self, chara):
        if chara == self.characters[-1]:
            return len(self.lc)
        else:
            return self.first[self.characters[self.characters.index(chara) + 1]]


# BWTstore("asdfqwer").find("dfqw")
# BWTstore("asdfqsdagfs").find("dfqsd")
print(BWTstore("asdfasdf").find("asdf"))

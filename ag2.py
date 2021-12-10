from random import seed
from typing import List
import numpy as np
import generator
import ag_tra


class BWTstore:
    def __init__(self, seq) -> None:
        seq_r = seq[::-1]
        fc = np.array(list("$" + seq_r))
        lc = np.array(list(seq_r + "$"))

        characters = np.unique(fc)
        characters.sort()
        self.characters = characters
        group0 = {key: np.where(fc == key)[0] for key in characters}

        def my_sort(s, group=None, all_position=None):

            # return local transmation if recursion
            # if len(s)=len(fl) -> local=global

            if len(s) == 0:
                pass
            if len(s) == 1:
                return [0]

            if len(s) == len(np.unique(s)):
                if group:
                    return list(np.argsort(s))
                else:
                    return list(all_position[np.argsort(s)])

            transformation = []

            if not group:
                letters = np.unique(s)
                letters.sort()
            else:
                letters = self.characters
            for letter in letters:
                if group:
                    # fl, has group=group0
                    position = group[letter]
                else:
                    # subseq, has all_position
                    position = np.array(all_position)[np.where(s == letter)]

                if len(position) < 2:
                    transformation += list(position)
                else:
                    pre_position = position - 1
                    pre_s = fc[pre_position]
                    local_transformation = my_sort(pre_s, all_position=pre_position)

                    local_transformation = list(np.array(local_transformation) + 1)
                    transformation += local_transformation

            return transformation

        transformation = my_sort(fc, group=group0)

        self.po = len(fc) - np.array(transformation)
        self.lc = lc[transformation]

        long = len(fc)
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
        for i in self.characters:
            first[i] = start
            start += counter[i]

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
            return self.first[self.characters[list(self.characters).index(chara) + 1]]

    def seeding(self, motif, k):
        max_match = len(motif) // (k + 1)
        overlap = min(max(2, int(max_match * 0.2)), max_match - 1)
        seed_search = []
        for i in range((len(motif) - overlap) // (max_match - overlap)):

            if i == 0:
                seed = motif[0:max_match]
            else:
                seed = motif[
                    i * (max_match - overlap) : i * (max_match - overlap) + max_match
                ]
            seed_search += self.find(seed)
        return seed_search, max_match, overlap


# print(BWTstore("asdfasdf").find("asdf"))

if __name__ == "__main__":
    print(
        BWTstore(generator.generate_DNA(100000000)).seeding(
            generator.generate_DNA(80), 10
        )
    )
    # print(
    #     BWTstore(
    #         "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAASSSSSSSSSSSAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    #     ).seeding("SSSSS", 1)
    # )

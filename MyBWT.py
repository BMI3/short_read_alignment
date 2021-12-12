import numpy as np


class MyBWT:
    def __init__(self, ref) -> None:

        self.characters, self.po, self.lc = self.generateLC(ref)
        self.tally = self.generateTally()
        self.char_range = self.charRange()

    def generateLC(self, ref):

        ref_r = ref[::-1]
        fc = np.array(list("$" + ref_r))
        lc = np.array(list(ref_r + "$"))

        characters = sorted(np.unique(fc))

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
                letters = sorted(np.unique(s))
            else:
                letters = characters
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

        po = len(fc) - np.array(transformation)
        lc = lc[transformation]
        return characters, po, lc

    def generateTally(self):
        """
        Burrow-Wheeler transformed string, return a dict called all_char_tally,
        key: char, value: its count from start to end (a list)
        """
        tally = {}
        all_char_tally = {}
        for c in self.characters:
            tally[c] = 0
            all_char_tally[c] = []
        for char in self.lc:
            tally[char] += 1
            for c in self.characters:
                all_char_tally[c].append(tally[c])
        return all_char_tally

    def charRange(self):
        """Given total: a dict with chars(nucleotide) as keys and its total occurences as values,
        return the postion of each char in the F column (i.e. sorted bwt) in which range it occurs (a range
        with right closed, left open)"""

        total = {key: value[-1] for key, value in self.tally.items()}
        char_range = dict()
        pos = 0
        for char, count in sorted(total.items()):
            char_range[char] = (pos, pos + count)
            pos += count
        return char_range

    def query(self, shortRead):
        """
        Given sequence to be queried and origin string t, performing the BWT,
        return the position that seq occur in the origin string (-1 for not found)
        """
        # get the F col
        fc = sorted(self.lc)
        # get the start and end pos of each char
        char_range = self.char_range

        # start from the last char in seq
        first_char = shortRead[-1]
        # check whether char in the ref string first
        if first_char not in self.characters:
            return -1
        # the range of the first char can be found in char_range
        cc_range = [char_range[first_char][0], char_range[first_char][1]]
        # length of seq is 1 (uncommon case)
        if len(shortRead) == 1:
            return self.po[cc_range[0] : cc_range[1]]
        # for the remaining char in seq
        for i in range(len(shortRead) - 2, -1, -1):
            cc_range = self.findNextWithTally(cc_range, shortRead[i], char_range)
            # print(seq[i],cc_range)
            if cc_range == -1:
                return -1
        return self.po[cc_range[0] : cc_range[1]]

    def findNextWithTally(self, current_range, next_char, char_range):
        """
        Given the current range (left closed, right open) in the F col and the next char searching for, with the help of tally and char_range,
        return the rank of the next char satified (-1 for not found)
        """
        # get the tally for all chars
        tally = self.tally

        start = current_range[0] - 1
        # need to see on before so that we do not leave out the first char
        end = current_range[1] - 1
        # right open

        if next_char not in tally:
            return -1
        num_char_found = tally[next_char][end] - tally[next_char][start]
        # how many next char found
        if num_char_found > 0:
            rank_start, rank_end = tally[next_char][start], tally[next_char][end] - 1
            return [
                char_range[next_char][0] + rank_start,
                char_range[next_char][0] + rank_end + 1,
            ]  # range is left open, right closed
        else:
            return -1

    def seeding(self, shortRead, k):
        max_match = len(shortRead) // (k + 1)
        overlap = min(max(2, int(max_match * 0.2)), max_match - 1)
        seed_search = []
        for i in range((len(shortRead) - overlap) // (max_match - overlap)):

            if i == 0:
                seed = shortRead[0:max_match]
            else:
                seed = shortRead[
                    i * (max_match - overlap) : i * (max_match - overlap) + max_match
                ]
            seed_search += self.query(seed)
        return seed_search, max_match, overlap

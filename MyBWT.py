import numpy as np


class MyBWT:
    def __init__(self, id, ref) -> None:
        self.id = id
        self.characters, self.po, self.lc = self.generateLC(ref)
        self.tally = self.generateTally()
        self.char_range = self.charRange()
        self.origin = ref

    def generateLC(self, ref):
        """
        generate the last characters of rotation in BWT and also record the positions in the origin string
        """
        ref_r = ref[::-1]
        fc = np.array(list("$" + ref_r))
        lc = np.array(list(ref_r + "$"))

        characters = sorted(np.unique(fc))

        group0 = {key: np.where(fc == key)[0] for key in characters}

        def my_sort(s, group=None, all_position=None):
            """
            return local transformation during recursion
            if s=fc -> local=global
            when s=fc -> group=True
                 len(s)<fc (deal subgroups) ->use all position
            """
            if len(s) == 0:
                pass
            if len(s) == 1:
                return [0]

            # if the string dose not contain repeat character, just argsort()
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

                if len(position) < 2:  # must be unique
                    transformation += list(position)
                else:
                    pre_position = position - 1
                    pre_s = fc[pre_position]  # the previous character in fc
                    local_transformation = my_sort(pre_s, all_position=pre_position)
                    local_transformation = list(np.array(local_transformation) + 1)
                    transformation += local_transformation

            return transformation

        transformation = my_sort(fc, group=group0)

        po = len(fc) - np.array(transformation)  # position in the origin string
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
        """
        Given total: a dict with chars(nucleotide) as keys and its total occurences as values,
        return the postion of each char in the F column (i.e. sorted bwt) in which range it occurs (a range
        with right closed, left open)
        """
        char_range = {}
        pos = 0
        for char, list in sorted(self.tally.items()):
            char_range[char] = (pos, pos + list[-1])
            pos += list[-1]
        return char_range

    def query(self, shortRead):
        """
        Given sequence to be queried and origin string t, performing the BWT,
        return the position that seq occur in the origin string (-1 for not found)
        """
        # get the start and end pos of each char
        char_range = self.char_range

        # start from the last char in seq
        first_char = shortRead[-1]
        # check whether char in the ref string first
        if first_char not in self.characters:
            return np.array([])
        # the range of the first char can be found in char_range
        cc_range = [char_range[first_char][0], char_range[first_char][1]]

        def findNextWithTally(current_range, next_char, char_range):
            """
            Given the current range (left closed, right open) in the F col and the next char searching for, with the help of tally and char_range,
            return the rank of the next char satified (-1 for not found)
            """
            tally = self.tally
            # get the tally for all chars
            start = current_range[0] - 1
            # need to see on before so that we do not leave out the first char
            end = min(current_range[1] - 1, self.lc.size - 2)
            # right open

            if next_char not in tally:
                return []
            num_char_found = tally[next_char][end] - tally[next_char][start]
            # how many next char found
            if num_char_found > 0:
                # rank_start, rank_end = tally[next_char][start], tally[next_char][end] - 1
                rank_start, rank_end = tally[next_char][start], tally[next_char][end]
                return (
                    char_range[next_char][0] + rank_start,
                    char_range[next_char][0] + rank_end + 1,
                )  # range is left open, right closed
            else:
                return []

        # for the remaining char in seq
        for i in range(len(shortRead) - 2, -1, -1):
            cc_range = findNextWithTally(cc_range, shortRead[i], char_range)
            if cc_range == []:
                return np.array([])
        return self.po[cc_range[0] : cc_range[1]]

    def seeding(self, shortRead, k=2):
        """
        split short reads into pieces and do query
        query will return the position of seeds
        convert the possition by minus the distance of the start of the seed in the short read
        therefore, all the seeds should have the same return position if no varience, insert, delete
        """
        max_match = len(shortRead) // (k + 1)
        overlap = min(max(2, int(max_match * 0.2)), max_match - 1)
        seed_search = []
        n_seed = (len(shortRead) - overlap) // (max_match - overlap) + 1
        for i in range(n_seed):
            seed = shortRead[
                i
                * (max_match - overlap) : min(
                    i * (max_match - overlap) + max_match, len(shortRead)
                )
            ]
            result = self.query(seed)
            if result.size > 0:
                start = result - i * (max_match - overlap)
                if len(np.where(start < 0)[0]):
                    start = start[np.where(start >= 0)]
                seed_search += list(start)

        if not seed_search:
            return seed_search

        possible = []
        freq = []
        for p in seed_search:
            if p not in possible:
                possible.append(p)
                freq.append(1)
            else:
                freq[possible.index(p)] += 1
        # print("possible", possible)
        # print("freq:", freq)
        return self.extend(shortRead, possible)
        # if max(freq) == n_seed:
        #     position = [possible[i] for i in range(len(freq)) if freq[i] == max(freq)]
        #     # if 680580 in position:
        #     #     raise error

        #     # self.extend(shortRead, position)

        #     return position
        # else:
        #     return possible

    # def traceback(self):
    #     """
    #     Given Burrows-Wheeler tranformed string bwt, return
    #     the origin sequence
    #     """
    #     rownum = 0
    #     seq = self.lc[0]

    #     for i in range(self.lc.size - 1):
    #         letter = self.lc[rownum]  # record the focus character in the last column
    #         rank = self.tally[letter][rownum]
    #         rownum = self.char_range[letter][0] + rank - 1
    #         seq = letter + seq
    #     self.origin = seq

    def extend(self, shortRead, position, threshold=0.95):
        """
        simply compare the part of origin reference on the returned position
                    with the shortread
        use Hamming distance to get a match score
        return position with score higher than the threshold
        """
        acceptable = []
        for p in position:
            p -= 1
            # if self.origin == None:
            #     self.traceback()

            o = self.origin[p : p + len(shortRead)]
            if len(o) != len(shortRead):
                sl = min(len(o), len(shortRead))
                o = o[:sl]
                shortRead = shortRead[:sl]
            match = sum((np.array(list(o)) == np.array(list(shortRead)))) / len(
                shortRead
            )
            if match > threshold:
                acceptable.append((p, match))
                # print(p, ":", match / len(shortRead) * 100, "%")
        return acceptable

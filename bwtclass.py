import numpy as np

class BWTtrans:

    '''Manage Burrow-Wheeler transform, query'''

    def __init__(self, ref, seq) -> None:
        '''create bwt string'''
        self.ref = ref
        self.seq = seq
        self.sa = self.suffixArray(self.ref)
        self.bwt = self.bwaBySa(self.ref,self.sa)
        self.total, self.ranks = self.trank()
        self.char_range = self.charRange(self.total)
        self.all_tally = self.generateTally()


    def suffixArray(self, t):
        ''' given string t, return its suffix array'''
        t = t + "$"
        sa = []
        for i in range(len(t)):
            sa.append((t[i:],i))
        return sorted(sa)

    def bwaBySa(self, t,sa = None):
        ''' given string t, return Burrow-Wheeler transformed string from the suffix array'''
        bwt = ''
        if sa == None:
            sa = self.suffixArray(t)
        for i in sa:
            if i[1] == 0:
                bwt += "$"
            else:
                bwt += t[i[1]-1]
        return bwt

    def trank(self):
        '''Given Burrows-Wheeler tranformed string, return
        total: a dict with chars(nucleotide) as keys and its total occurences as values
        ranks: a list of the char's rank'''
        bwt = self.bwt
        total = dict()
        ranks = []
        for i in bwt:
            if i not in total:
                total[i] = 0
            ranks.append(total[i])
            total[i] += 1
        return total, ranks

    def charRange(total):
        '''Given total: a dict with chars(nucleotide) as keys and its total occurences as values,
        return the postion of each char in the F column (i.e. sorted bwt) in which range it occurs (a range
        with right closed, left open)'''
        char_range = dict()
        pos = 0
        for char, count in sorted(total.items()):
            char_range[char] = (pos,pos + count)
            pos += count
        return char_range

    def findFPos(Lpos,char,ranks,start_pos):
        '''Given the pos in the L col (Lpos) and char, with ranks from trank() and start_pos
        from startPos(), return a list next_pos which is the pos in the F col that matches
        the char in the L col(like traceback)'''
        next_pos = []
        for p in Lpos:
            # ranks[p]: the rank we can jump; start_pos: the char we can jump
            next_pos.append(ranks[p]+start_pos[char])
        return next_pos

    def generateTally():
        ''' Burrow-Wheeler transformed string, return a dict called all_char_tally,
        key: char, value: its count from start to end (a list)'''
        bwt = self.bwt
        tally = dict()
        all_char_tally = dict()
        for i in bwt:
            if i not in tally:
                tally[i] = 0
                all_char_tally[i] = []
        for rowindex, char in enumerate(bwt):
            tally[char] += 1
            for c in tally.keys():
                all_char_tally[c].append(tally[c])
        return all_char_tally

    def findNextWithTally(current_range, next_char,tally,char_range):
        '''Given the current range (left closed, right open) in the F col and the next char searching for, with the help of tally and char_range,
        return the rank of the next char satified (-1 for not found)'''
        start = current_range[0] - 1 # need to see on before so that we do not leave out the first char
        end = current_range[1] - 1 # right open
        if next_char not in tally:
            return -1
        num_char_found = tally[next_char][end] - tally[next_char][start] # how many next char found
        if num_char_found > 0:
            rank_start, rank_end = tally[next_char][start], tally[next_char][end]-1
            return [char_range[next_char][0]+rank_start, char_range[next_char][0]+rank_end+1] # range is left open, right closed
        else:
            return -1
    
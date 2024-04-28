from Bio.Align import substitution_matrices
import numpy as np

def Needleman_Wunsch(seq1, seq2, del_cost, ins_cost, mismatch_penalty, match_score, matrix = False):

    m, n = len(seq1), len(seq2)
    dp = np.zeros((m + 1, n + 1))
    dp[:, 0] = del_cost * np.arange(m + 1)
    dp[0, :] = ins_cost * np.arange(n + 1)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == '-' and seq2[j - 1] != "-":
                score = del_cost
            elif seq2[j - 1] == '-' and seq1[i - 1] != "-":
                score = ins_cost
            elif seq1[i - 1] == seq2[j - 1] and seq1[i - 1] != "-":
                score = match_score
            else:
                score = mismatch_penalty
            dp[i, j] = max(
                dp[i - 1, j] + del_cost,
                dp[i, j - 1] + ins_cost,
                dp[i - 1, j - 1] + score,
            )
    if matrix == False:
        return dp[m][n]
    return dp

def get_alignment(seq1, seq2, del_cost, ins_cost, mismatch_penalty, match_score, dp):
    align1, align2 = [], []
    i, j = len(seq1), len(seq2)
    while i + j > 0:
        if i > 0 and dp[i, j] == dp[i - 1, j] + ins_cost:
            i -= 1
            align1.append(seq1[i])
            align2.append('-')
        elif j > 0 and dp[i, j] == dp[i, j - 1] + del_cost:
            j -= 1
            align1.append('-')
            align2.append(seq2[j])
        else:
            i -= 1
            j -= 1
            align1.append(seq1[i])
            align2.append(seq2[j])

    return (''.join(align1[::-1]), ''.join(align2[::-1]))


class Consensus: # Даже последовательности являются консенсусами
    def __init__(self, seq, p = None, c = 1, res = None):
        idxs = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4}
        self.seq = seq
        if p is None :
            p = np.zeros((5, len(seq)))
            for i in range(len(seq)):
                p[idxs[seq[i]]][i] = 1
        self.profile = p
        self.count = c
        if res is None:
            self.result = np.array(list(seq))
        else:
            self.result = res


def merge(con1, con2, del_cost, ins_cost, match_score, mismatch_penalty):
    matrix = Needleman_Wunsch(con1.seq, con2.seq, del_cost, ins_cost, match_score, mismatch_penalty, True)
    seq1, seq2 = get_alignment(con1.seq, con2.seq, del_cost, ins_cost, match_score, mismatch_penalty, matrix)



    l = max(len(seq1), len(seq2))
    p = np.zeros((5, l))
    gap = 0
    for i, val in enumerate(seq1):
        if (val == "-" and i == len(con1.seq)) or val != con1.seq[i-gap]:
            con1.profile = np.insert(con1.profile, i-gap, [0, 0, 0, 0, con1.count], axis= 1)
            if con1.result.ndim == 1 :
                con1.result = np.insert(con1.result, i-gap, ["-" for i in range(con1.count)])
            else:
                con1.result = np.insert(con1.result, i-gap, ["-" for i in range(con1.count)], axis= 1)
            gap+=1

    gap = 0
    for i, val in enumerate(seq2):
        if (val == "-" and i == len(con2.seq)) or val != con2.seq[i-gap]:
            con2.profile = np.insert (con2.profile, i-gap, [0, 0, 0, 0, con2.count], axis= 1)
            if con2.result.ndim == 1:
                con2.result = np.insert(con2.result, i-gap, ["-" for i in range(con2.count)])
            else:
                con2.result = np.insert(con2.result, i-gap, ["-" for i in range(con2.count)], axis= 1)
            gap+=1

    p = con1.profile + con2.profile
    c = con1.count + con2.count
    new_con_int = np.argmax(p, axis=0)
    revert_idxs = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: '-'}
    new_con = ""

    res_new = np.row_stack([con2.result, con1.result])
    for val in new_con_int:
        new_con += revert_idxs[val]
    new = Consensus(new_con, p, c, res_new)
    return new

def greedy_multiple_alignment(seq, del_cost, ins_cost, match_cost, mismatch_cost):
    l = len(seq)
    lst = []
    for i in range(l):
        lst.append(Consensus(seq[i].upper()))
    matrix = np.zeros((l, l))

    for i in range(l):
        for j in range(l):
            if i == j :
                matrix[i][j] = np.nan
            else:
                matrix[i][j] = Needleman_Wunsch(lst[i].seq, lst[j].seq, del_cost, ins_cost, match_cost, mismatch_cost)
    for k in range(l-1):
        ri, ci = np.nanargmax(matrix)//matrix.shape[1], np.nanargmax(matrix)%matrix.shape[1]
        lst[ri] = merge(lst[ri], lst[ci], del_cost, ins_cost, match_cost, mismatch_cost)
        lst[ci] = None
        matrix[:,ci] = np.nan
        matrix[ci,:] = np.nan
        for i in range(l):
            if not np.isnan(matrix[i,ri]):
                matrix[i,ri] = Needleman_Wunsch(lst[i].seq, lst[ri].seq, del_cost, ins_cost, match_cost, mismatch_cost)
            if not np.isnan(matrix[ri,i]):
                matrix[ri,i] = Needleman_Wunsch(lst[ri].seq, lst[i].seq, del_cost, ins_cost, match_cost, mismatch_cost)
    for i in range(l):
        if lst[i] != None:
            L = lst[i].result.tolist()
            for k in L:
                print(" ".join(k))
            return lst[i].seq

if __name__ == '__main__':
    seqs_1 = ["ACT", "ATC", "GCT", "ATCC"]
    seqs_2 = [
        'TTGGGGACTTCC',
        'TCGGGGATTCAT',
        'TCGGGGATTCCT',
        'TAGGGGAACTAC',
        'TCGGGTATAACC',
        'TCGGGGGTTTTT',
        'CCGGTGACTTAC',
        'ACGGGGATTTTC',
        'TTGGGGACTTTT',
        'AAGGGGACTTCC',
    ]
    seqs_3 = [
        "ttgacagctagctcagtcctaggtataatgctagc",
        "ttgacagctagctcagtcctaggtataatgctagc",
        "tttacagctagctcagtcctaggtattatgctagc",
        "ttgacagctagctcagtcctaggtactgtgctagc",
        "ctgatagctagctcagtcctagggattatgctagc",
        "ttgacagctagctcagtcctaggtattgtgctagc",
        "tttacggctagctcagtcctaggtactatgctagc",
        "tttacggctagctcagtcctaggtatagtgctagc",
        "tttacggctagctcagccctaggtattatgctagc",
        "ctgacagctagctcagtcctaggtataatgctagc",
        "tttacagctagctcagtcctagggactgtgctagc",
        "tttacggctagctcagtcctaggtacaatgctagc",
        "ttgacggctagctcagtcctaggtatagtgctagc",
        "ctgatagctagctcagtcctagggattatgctagc",
        "ctgatggctagctcagtcctagggattatgctagc",
        "tttatggctagctcagtcctaggtacaatgctagc",
        "tttatagctagctcagcccttggtacaatgctagc",
        "ttgacagctagctcagtcctagggactatgctagc",
        "ttgacagctagctcagtcctagggattgtgctagc",
        "ttgacggctagctcagtcctaggtattgtgctagc"
    ]
    ans = greedy_multiple_alignment(seqs_2,-2,-2,-1,1)
    print("Консенсусная последовательность", ans)


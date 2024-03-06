from Bio.Align import substitution_matrices

##################### Task 1 practice
def Needleman_Wunsch(seq1, seq2, matrix, gap):

    def score(pair):
        if pair in matrix:
            return matrix[pair]
        else:
            return matrix[tuple(reversed(pair))]

    dp = [[0 for i in range(len(seq1)+1)] for j in range(len(seq2)+1)]
    pointers = [["" for i in range(len(seq1)+1)] for j in range(len(seq2)+1)]

    dp[0] = [-gap*i for i in range(len(seq1)+1)]
    for j,row in enumerate(dp): #base dp
        row[0] = -gap*float(j)


    for i in range(len(seq2)+1)[1:]:
        for j in range(len(seq1)+1)[1:]:

            f = [dp[i-1][j-1] + score((seq2[i-1],seq1[j-1])), # match
                 dp[i-1][j] - gap,                              # gap
                 dp[i][j-1] - gap]                              # gap

            dp[i][j] = max(f)
            argmax = f.index(max(f))
            pointers[i][j] = argmax

    final_score = dp[len(seq2)][len(seq1)]

    align1 = ""
    align2 = ""
    i = len(seq2)
    j = len(seq1)

    while True:
        if i==0 or j==0:
            break


        if pointers[i][j] == 0: #diagonal
            align1 = seq1[j-1] + align1
            align2 = seq2[i-1] + align2

            i = i-1
            j = j-1
        elif pointers[i][j] == 1: #up
            align1 = "-" + align1
            align2 = seq2[i-1] + align2

            i = i-1

        elif pointers[i][j] == 2: #left
            align1 = seq1[j-1] + align1
            align2 = "-" + align2

            j = j-1

    print(align1)
    print(align2)
    print(final_score)



matrix = substitution_matrices.load('BLOSUM62')
seq1 = "GTTAC"
seq2 = "GACGT"
gap = 3
Needleman_Wunsch(seq1, seq2, matrix, gap)

##################### Task 2 practice
def Needleman_Wunsch_affinities(seq1, seq2, matrix, alpha, beta):

    def score(pair):
        if pair in matrix:
            return matrix[pair]
        else:
            return matrix[tuple(reversed(pair))]

    dp = [[0 for i in range(len(seq1)+1)] for j in range(len(seq2)+1)]
    pointers = [["" for i in range(len(seq1)+1)] for j in range(len(seq2)+1)]

    dp[0] = [-alpha*i for i in range(len(seq1)+1)]
    for j,row in enumerate(dp): #base dp
        row[0] = -alpha*float(j)


    for i in range(len(seq2)+1)[1:]:
        for j in range(len(seq1)+1)[1:]:
            gap_first = alpha
            gap_second = alpha
            if pointers[i-1][j] == 1:
                gap_first = beta
            if pointers[i][j-1] == 2:
                gap_second = beta

            f = [dp[i-1][j-1] + score((seq2[i-1],seq1[j-1])), # match
                 dp[i-1][j] - gap_first,                              # gap
                 dp[i][j-1] - gap_second]                              # gap

            dp[i][j] = max(f)
            argmax = f.index(max(f))
            pointers[i][j] = argmax

    final_score = dp[len(seq2)][len(seq1)]

    align1 = ""
    align2 = ""
    i = len(seq2)
    j = len(seq1)

    while True:
        if i==0 or j==0:
            break


        if pointers[i][j] == 0: #diagonal
            align1 = seq1[j-1] + align1
            align2 = seq2[i-1] + align2

            i = i-1
            j = j-1
        elif pointers[i][j] == 1: #up
            align1 = "-" + align1
            align2 = seq2[i-1] + align2

            i = i-1

        elif pointers[i][j] == 2: #left
            align1 = seq1[j-1] + align1
            align2 = "-" + align2

            j = j-1

    print(align1)
    print(align2)
    print(final_score)


matrix = substitution_matrices.load('BLOSUM62')
seq1 = "GTTAC"
seq2 = "GACGT"
alpha = 3
beta = 1
Needleman_Wunsch_affinities(seq1, seq2, matrix, alpha, beta)

##################### Theoretical part
""" 1. Количество выравниваний можно получить суммой выравниваний при выпадении последнего нуклеотида, либо вставке
    в зависимости от рассматриваемой строки. При этом ситуацию совпадения мы не рассматриваем, так как она может осуществиться либо
    в первый случай, либо во второй и не изменяет количество - тогда рекуррентная формула:
    N(n, m) = N(n - 1, m) + N(n, m - 1) 
    
    2. По определению количество выравниваний - формула сочетаний из комбинаторики. Докажем это по матиндукции:
    N(1, 0) = 1 
    N(0, 1) = 1 
    N(1, 1) = N(1, 0) + N(0, 1) = 2
    В предположении, что N(n - 1, m) = (n+m-1)!/((m)!(n-1)!), N(n, m - 1) = (n+m-1)!/((m-1)!n!
    N(n, m) = N(n - 1, m) + N(n, m - 1) = (n+m-1)!/((m)!(n-1)!) + (n+m-1)!/((m-1)!n! = (n+m)!/(m!n!), что
    является числом сочетаний из m+n по m. чтд
    
    3. Приближение Стирлинга :
    n! ~ sqrt(2 pi n) (n / e)^n Допустим, что n~=m, тогда: 
    N(n, m) = (n+m)!/(m!m!) = 1/sqrt(2pi) * sqrt(2pi) * (2n)^2n / (n^2n+1) = 2^(2n) / sqrt(pi*n)
    
    """




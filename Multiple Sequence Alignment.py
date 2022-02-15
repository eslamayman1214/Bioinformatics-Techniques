import numpy as nm
blosum62 = {
    '-':{'-':1,'A':-4,'C':-4,'B':-4,'E':-4,'D':-4,'G':-4,'F':-4,'I':-4,'H':-4,'K':-4,'M':-4,'L':-4,'N':-4,'Q':-4,'P':-4,'S':-4,'R':-4,'T':-4,'W':-4,'V':-4,'Y':-4,'X':-4,'Z':-4},
    'A':{'-':-4,'A':4,'C':0,'B':-2,'E':-1,'D':-2,'G':0,'F':-2,'I':-1,'H':-2,'K':-1,'M':-1,'L':-1,'N':-2,'Q':-1,'P':-1,'S':1,'R':-1,'T':0,'W':-3,'V':0,'Y':-2,'X':-1,'Z':-1},
    'C':{'-':-4,'A':0,'C':9,'B':-3,'E':-4,'D':-3,'G':-3,'F':-2,'I':-1,'H':-3,'K':-3,'M':-1,'L':-1,'N':-3,'Q':-3,'P':-3,'S':-1,'R':-3,'T':-1,'W':-2,'V':-1,'Y':-2,'X':-1,'Z':-3},
    'B':{'-':-4,'A':-2,'C':-3,'B':4,'E':1,'D':4,'G':-1,'F':-3,'I':-3,'H':0,'K':0,'M':-3,'L':-4,'N':3,'Q':0,'P':-2,'S':0,'R':-1,'T':-1,'W':-4,'V':-3,'Y':-3,'X':-1,'Z':1},
    'E':{'-':-4,'A':-1,'C':-4,'B':1,'E':5,'D':2,'G':-2,'F':-3,'I':-3,'H':0,'K':1,'M':-2,'L':-3,'N':0,'Q':2,'P':-1,'S':0,'R':0,'T':-1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':4},
    'D':{'-':-4,'A':-2,'C':-3,'B':4,'E':2,'D':6,'G':-1,'F':-3,'I':-3,'H':-1,'K':-1,'M':-3,'L':-4,'N':1,'Q':0,'P':-1,'S':0,'R':-2,'T':-1,'W':-4,'V':-3,'Y':-3,'X':-1,'Z':1},
    'G':{'-':-4,'A':0,'C':-3,'B':-1,'E':-2,'D':-1,'G':6,'F':-3,'I':-4,'H':-2,'K':-2,'M':-3,'L':-4,'N':0,'Q':-2,'P':-2,'S':0,'R':-2,'T':-2,'W':-2,'V':-3,'Y':-3,'X':-1,'Z':-2},
    'F':{'-':-4,'A':-2,'C':-2,'B':-3,'E':-3,'D':-3,'G':-3,'F':6,'I':0,'H':-1,'K':-3,'M':0,'L':0,'N':-3,'Q':-3,'P':-4,'S':-2,'R':-3,'T':-2,'W':1,'V':-1,'Y':3,'X':-1,'Z':-3},
    'I':{'-':-4,'A':-1,'C':-1,'B':-3,'E':-3,'D':-3,'G':-4,'F':0,'I':4,'H':-3,'K':-3,'M':1,'L':2,'N':-3,'Q':-3,'P':-3,'S':-2,'R':-3,'T':-1,'W':-3,'V':3,'Y':-1,'X':-1,'Z':-3},
    'H':{'-':-4,'A':-2,'C':-3,'B':0,'E':0,'D':-1,'G':-2,'F':-1,'I':-3,'H':8,'K':-1,'M':-2,'L':-3,'N':1,'Q':0,'P':-2,'S':-1,'R':0,'T':-2,'W':-2,'V':-3,'Y':2,'X':-1,'Z':0},
    'K':{'-':-4,'A':-1,'C':-3,'B':0,'E':1,'D':-1,'G':-2,'F':-3,'I':-3,'H':-1,'K':5,'M':-1,'L':-2,'N':0,'Q':1,'P':-1,'S':0,'R':2,'T':-1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':1},
    'M':{'-':-4,'A':-1,'C':-1,'B':-3,'E':-2,'D':-3,'G':-3,'F':0,'I':1,'H':-2,'K':-1,'M':5,'L':2,'N':-2,'Q':0,'P':-2,'S':-1,'R':-1,'T':-1,'W':-1,'V':1,'Y':-1,'X':-1,'Z':-1},
    'L':{'-':-4,'A':-1,'C':-1,'B':-4,'E':-3,'D':-4,'G':-4,'F':0,'I':2,'H':-3,'K':-2,'M':2,'L':4,'N':-3,'Q':-2,'P':-3,'S':-2,'R':-2,'T':-1,'W':-2,'V':1,'Y':-1,'X':-1,'Z':-3},
    'N':{'-':-4,'A':-2,'C':-3,'B':3,'E':0,'D':1,'G':0,'F':-3,'I':-3,'H':1,'K':0,'M':-2,'L':-3,'N':6,'Q':0,'P':-2,'S':1,'R':0,'T':0,'W':-4,'V':-3,'Y':-2,'X':-1,'Z':0},
    'Q':{'-':-4,'A':-1,'C':-3,'B':0,'E':2,'D':0,'G':-2,'F':-3,'I':-3,'H':0,'K':1,'M':0,'L':-2,'N':0,'Q':5,'P':-1,'S':0,'R':1,'T':-1,'W':-2,'V':-2,'Y':-1,'X':-1,'Z':3},
    'P':{'-':-4,'A':-1,'C':-3,'B':-2,'E':-1,'D':-1,'G':-2,'F':-4,'I':-3,'H':-2,'K':-1,'M':-2,'L':-3,'N':-2,'Q':-1,'P':7,'S':-1,'R':-2,'T':-1,'W':-4,'V':-2,'Y':-3,'X':-1,'Z':-1},
    'S':{'-':-4,'A':1,'C':-1,'B':0,'E':0,'D':0,'G':0,'F':-2,'I':-2,'H':-1,'K':0,'M':-1,'L':-2,'N':1,'Q':0,'P':-1,'S':4,'R':-1,'T':1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':0},
    'R':{'-':-4,'A':-1,'C':-3,'B':-1,'E':0,'D':-2,'G':-2,'F':-3,'I':-3,'H':0,'K':2,'M':-1,'L':-2,'N':0,'Q':1,'P':-2,'S':-1,'R':5,'T':-1,'W':-3,'V':-3,'Y':-2,'X':-1,'Z':0},
    'T':{'-':-4,'A':0,'C':-1,'B':-1,'E':-1,'D':-1,'G':-2,'F':-2,'I':-1,'H':-2,'K':-1,'M':-1,'L':-1,'N':0,'Q':-1,'P':-1,'S':1,'R':-1,'T':5,'W':-2,'V':0,'Y':-2,'X':-1,'Z':-1},
    'W':{'-':-4,'A':-3,'C':-2,'B':-4,'E':-3,'D':-4,'G':-2,'F':1,'I':-3,'H':-2,'K':-3,'M':-1,'L':-2,'N':-4,'Q':-2,'P':-4,'S':-3,'R':-3,'T':-2,'W':11,'V':-3,'Y':2,'X':-1,'Z':-3},
    'V':{'-':-4,'A':0,'C':-1,'B':-3,'E':-2,'D':-3,'G':-3,'F':-1,'I':3,'H':-3,'K':-2,'M':1,'L':1,'N':-3,'Q':-2,'P':-2,'S':-2,'R':-3,'T':0,'W':-3,'V':4,'Y':-1,'X':-1,'Z':-2},
    'Y':{'-':-4,'A':-2,'C':-2,'B':-3,'E':-2,'D':-3,'G':-3,'F':3,'I':-1,'H':2,'K':-2,'M':-1,'L':-1,'N':-2,'Q':-1,'P':-3,'S':-2,'R':-2,'T':-2,'W':2,'V':-1,'Y':7,'X':-1,'Z':-2},
    'X':{'-':-4,'A':-1,'C':-1,'B':-1,'E':-1,'D':-1,'G':-1,'F':-1,'I':-1,'H':-1,'K':-1,'M':-1,'L':-1,'N':-1,'Q':-1,'P':-1,'S':-1,'R':-1,'T':-1,'W':-1,'V':-1,'Y':-1,'X':-1,'Z':-1},
    'Z':{'-':-4,'A':-1,'C':-3,'B':1,'E':4,'D':1,'G':-2,'F':-3,'I':-3,'H':0,'K':1,'M':-1,'L':-3,'N':0,'Q':3,'P':-1,'S':0,'R':0,'T':-1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':4}}


def Unite_length_of_sequences(first_seq, secod_seq, third_seq, fourth_seq):
    if len(first_seq) < len(secod_seq):
        for i in range(len(first_seq) + 1, len(secod_seq) + 1):
            first_seq = first_seq + "-"

    if len(secod_seq) < len(first_seq):
        for i in range(len(secod_seq) + 1, len(first_seq) + 1):
            secod_seq = secod_seq + "-"

    if len(third_seq) < len(fourth_seq):
        for i in range(len(third_seq) + 1, len(fourth_seq) + 1):
            third_seq = third_seq + "-"

    if len(fourth_seq) < len(third_seq):
        for i in range(len(fourth_seq) + 1, len(third_seq) + 1):
            fourth_seq = fourth_seq + "-"

    return first_seq, secod_seq, third_seq, fourth_seq


match=1
missmatch=-1
gap=-1
def getmsascore(firstseq,secondseq,thirdseq,fourthseq):
    listofscores=[]
    score=0
    sequence1 =firstseq.upper()
    sequence2 =secondseq.upper()
    sequence3 =thirdseq.upper()
    sequence4 =fourthseq.upper()
    for i in range(0,len(sequence1)):
        for j in range(0,len(sequence3)):
            if sequence1[i]==sequence2[i]:
                if sequence1[i] == "-":
                    # handle gap gap score=0
                    score += 0
                else:
                    score+=match
                if sequence1[i]==sequence3[j]:
                    if sequence1[i] == "-":
                        # handle gap gap score=0
                        score += 0
                    else:
                        score += match
                    if sequence1[i]==sequence4[j]:
                        if sequence1[i] == "-":
                            # handle gap gap score=0
                            score += 0
                        else:
                            score += match
                    elif sequence1[i]=="-" or sequence4[j]=="-":
                        score+=gap
                    else :
                        score+=missmatch
                elif sequence1[i]=="-" or sequence3[j]=="-":
                    score+=gap
                    if sequence1[i]==sequence4[j]:
                        if sequence1[i] == "-":
                            # handle gap gap score=0
                            score += 0
                        else:
                            score += match
                    elif sequence1[i]=="-" or sequence4[j]=="-":
                        score+=gap
                    else :
                        score+=missmatch
                else:
                    score+=missmatch
                    if sequence1[i]==sequence4[j]:
                        if sequence1[i] == "-":
                            # handle gap gap score=0
                            score += 0
                        else:
                            score += match
                    elif sequence1[i]=="-" or sequence4[j]=="-":
                        score+=gap
                    else :
                        score+=missmatch
            elif sequence1[i]=="-" or sequence2[i]=="-":
                score+=gap
                if sequence1[i]==sequence3[j]:
                    if sequence1[i] == "-":
                        # handle gap gap score=0
                        score += 0
                    else:
                        score += match
                    if sequence1[i]==sequence4[j]:
                        if sequence1[i] == "-":
                            # handle gap gap score=0
                            score += 0
                        else:
                            score += match
                    elif sequence1[i]=="-" or sequence4[j]=="-":
                        score+=gap
                    else :
                        score+=missmatch
                elif sequence1[i]=="-" or sequence3[j]=="-":
                    score+=gap
                    if sequence1[i]==sequence4[j]:
                        if sequence1[i] == "-":
                            # handle gap gap score=0
                            score += 0
                        else:
                            score += match
                    elif sequence1[i]=="-" or sequence4[j]=="-":
                        score+=gap
                    else :
                        score+=missmatch
                else:
                    score+=missmatch
                    if sequence1[i]==sequence4[j]:
                        if sequence1[i] == "-":
                            # handle gap gap score=0
                            score += 0
                        else:
                            score += match
                    elif sequence1[i]=="-" or sequence4[j]=="-":
                        score+=gap
                    else :
                        score+=missmatch
            else:
                score+=missmatch
                if sequence1[i]==sequence3[j]:
                    if sequence1[i] == "-":
                        # handle gap gap score=0
                        score += 0
                    else:
                        score += match
                    if sequence1[i]==sequence4[j]:
                        if sequence1[i] == "-":
                            # handle gap gap score=0
                            score += 0
                        else:
                            score += match
                    elif sequence1[i]=="-" or sequence4[j]=="-":
                        score+=gap
                    else :
                        score+=missmatch
                elif sequence1[i]=="-" or sequence3[j]=="-":
                    score+=gap
                    if sequence1[i]==sequence4[j]:
                        if sequence1[i] == "-":
                            # handle gap gap score=0
                            score += 0
                        else:
                            score += match
                    elif sequence1[i]=="-" or sequence4[j]=="-":
                        score+=gap
                    else :
                        score+=missmatch
                else:
                    score+=missmatch
                    if sequence1[i]==sequence4[j]:
                        if sequence1[i] == "-":
                            # handle gap gap score=0
                            score += 0
                        else:
                            score += match
                    elif sequence1[i]=="-" or sequence4[j]=="-":
                        score+=gap
                    else :
                        score+=missmatch
            #second compare
            if sequence2[i]==sequence3[j]:
                if sequence2[i] == "-":
                    # handle gap gap score=0
                    score += 0
                else:
                    score += match
                if sequence2[i]==sequence4[j]:
                    if sequence2[i] == "-":
                        # handle gap gap score=0
                        score += 0
                    else:
                        score += match
                elif sequence2[i]=="-" or sequence4[j]:
                    score+=gap
                else:
                    score+=missmatch
            elif sequence2[i]=="-" or sequence3[j]=="-":
                score+=gap
                if sequence2[i]==sequence4[j]:
                    if sequence2[i] == "-":
                        # handle gap gap score=0
                        score += 0
                    else:
                        score += match
                elif sequence2[i]=="-" or sequence4[j]:
                    score+=gap
                else:
                    score+=missmatch
            else:
                score+=missmatch
                if sequence2[i]==sequence4[j]:
                    if sequence2[i] == "-":
                        # handle gap gap score=0
                        score += 0
                    else:
                        score += match
                elif sequence2[i]=="-" or sequence4[j]:
                    score+=gap
                else:
                    score+=missmatch
            # third compare
            if sequence3[j]==sequence4[j]:
                if sequence3[j] == "-":
                    # handle gap gap score=0
                    score += 0
                else:
                    score += match
            elif sequence3[j]=="-" or sequence4[j]=="-" :
                score+=gap
            else:
                score+=missmatch
            listofscores.append(score)
            score=0

    return listofscores




#firstseq is the orignal seq to call get score function
#sequence1 is " "+orignal sequence to implemen the matrix

def MSA_DNA(fisrtseq,secondseq,thirdseq,fourthseq,GAP):
    fisrtseq,secondseq,thirdseq,fourthseq=Unite_length_of_sequences(fisrtseq,secondseq,thirdseq,fourthseq)
    sequence1 = " " + fisrtseq.upper()
    sequence2 = " " + secondseq.upper()
    sequence3 = " " + thirdseq.upper()
    sequence4 = " " + fourthseq.upper()
    FirstAlign = []
    SecondAlign = []
    ThirdAlign = []
    FourthAlign = []
    NextStep = ""
    # k and n is the length of row and colum
    #to tracebak the alignment
    # k for sequence 1 and 2
    #n for sequence 3 and 4
    k =len(sequence1)-1
    n =len(sequence3)-1
    Gappenality=int(GAP)
    # initialize first row and colum
    Matrix = nm.zeros((len(sequence1), len(sequence3)))
    for i in range(len(sequence1)):
        Matrix[i][0] = i * Gappenality
    for j in range(len(sequence3)):
        Matrix[0][j] = j * Gappenality
# this list carry all scores of matrix before calulat left and up score diagonul score
# without add this score to last i-1 and j-1 score in matrix
    list_of_score=getmsascore(fisrtseq,secondseq,thirdseq,fourthseq)
    index_of_list=0
    for i in range(1, len(sequence1)):
         for j in range(1, len(sequence3)):
                 D = list_of_score[index_of_list] + Matrix[i - 1][j - 1]
                 L = Gappenality + Matrix[i][j - 1]
                 U = Gappenality + Matrix[i - 1][j]
                 MaxScoreIndex = nm.argmax([D, L, U])    # Return First Max index
                 if (MaxScoreIndex == 0):
                     NextStep += "D"
                 elif (MaxScoreIndex == 1):
                     NextStep += "L"
                 else:
                     NextStep += "U"
                 index_of_list += 1
                 Matrix[i][j] = nm.max([D, L, U])

    # Matrix Traceback
    Backward_Directions = nm.reshape(list(NextStep), (len(sequence1)-1 , len(sequence3)-1 ))
    Backward_Directions = nm.vstack([["*"] * Backward_Directions.shape[1], Backward_Directions])
    Backward_Directions = nm.column_stack([["*"] * Backward_Directions.shape[0], Backward_Directions])

    while (True):
        if (Backward_Directions[k][n] == "D"):
            FirstAlign.append(sequence1[k])
            SecondAlign.append(sequence2[k])
            ThirdAlign.append(sequence3[n])
            FourthAlign.append(sequence4[n])
            k -= 1
            n -= 1
        elif (Backward_Directions[k][n] == "U"):
            if (sequence1[k] == ' ' and sequence2[k]==' '):
                FirstAlign.append("-")
                SecondAlign.append("-")
                ThirdAlign.append(sequence3[n])
                FourthAlign.append(sequence4[n])
                n -= 1
            else:
                FirstAlign.append(sequence1[k])
                SecondAlign.append(sequence2[k])
                ThirdAlign.append("-")
                FourthAlign.append("-")
                k -= 1
        else:
            if (sequence3[n] == ' ' and sequence4[n]==' ' ):
                ThirdAlign.append("-")
                FourthAlign.append("-")
                FirstAlign.append(sequence1[k])
                SecondAlign.append(sequence2[k])
                k -= 1
            else:
                ThirdAlign.append(sequence3[n])
                FourthAlign.append(sequence4[n])
                FirstAlign.append("-")
                SecondAlign.append("-")
                n -= 1

        if (k < 0):
            break
        if (n < 0):
            break
    FirstAlign.reverse()
    SecondAlign.reverse()
    ThirdAlign.reverse()
    FourthAlign.reverse()

    print("The Best Alignment Sequences are :")
    print("Alignment Sequence 1 : " , FirstAlign[1:])
    print("Alignment Sequence 2 : " , SecondAlign[1:])
    print("Alignment Sequence 3 : ", ThirdAlign[1:])
    print("Alignment Sequence 4 : ", FourthAlign[1:])
    print("\nThe Scoring Matrix : \n")
    print(Matrix)
    print("The Score Of The Cell We begin From is  : ",  Matrix[len(sequence1) - 1][len(sequence3) - 1], "\n")



def getmsascore_for_protein(firstseq,secondseq,thirdseq,fourthseq):
    listofscores=[]
    score=0
    sequence1 =firstseq.upper()
    sequence2 =secondseq.upper()
    sequence3 =thirdseq.upper()
    sequence4 =fourthseq.upper()
    for i in range(0,len(sequence1)):
        for j in range(0,len(sequence3)):
            score += blosum62[sequence1[i]][sequence2[i]]
            score += blosum62[sequence1[i]][sequence3[j]]
            score += blosum62[sequence1[i]][sequence4[j]]
            score += blosum62[sequence2[i]][sequence3[j]]
            score += blosum62[sequence2[i]][sequence4[j]]
            score += blosum62[sequence3[j]][sequence4[j]]
            listofscores.append(score)
            score = 0

    return listofscores



def MSA_protein(fisrtseq,secondseq,thirdseq,fourthseq,gap):
    fisrtseq, secondseq, thirdseq, fourthseq = Unite_length_of_sequences(fisrtseq, secondseq, thirdseq, fourthseq)
    sequence1 = " " + fisrtseq.upper()
    sequence2 = " " + secondseq.upper()
    sequence3 = " " + thirdseq.upper()
    sequence4 = " " + fourthseq.upper()
    MScore = 1
    MMScore = -1
    GScore = -1
    FirstAlign = []
    SecondAlign = []
    ThirddAlign = []
    FourthAlign = []
    NextStep = ""
    k = len(sequence1) - 1
    n = len(sequence3) - 1
    GAP=int(gap)
    # initialize first row and colum
    Matrix = nm.zeros((len(sequence1), len(sequence3)))
    for i in range(len(sequence1)):
        Matrix[i][0] = i * GAP
    for j in range(len(sequence3)):
        Matrix[0][j] = j * GAP
    index_of_matrix = 0
    list_of_score = getmsascore_for_protein(fisrtseq,secondseq,thirdseq,fourthseq)

    #MATRIX CALCULATION
    for i in range(1, len(sequence1)):
        for j in range(1, len(sequence3)):
            D = list_of_score[index_of_matrix] + Matrix[i - 1][j - 1]
            L = GAP + Matrix[i][j - 1]
            U = GAP + Matrix[i - 1][j]
            MaxScoreIndex = nm.argmax([D, L, U])
            if (MaxScoreIndex == 0):
                NextStep += "D"
            elif (MaxScoreIndex == 1):
                NextStep += "L"
            else:
                NextStep += "U"
            index_of_matrix += 1
            Matrix[i][j] = nm.max([D, L, U])

    Backward_Directions = nm.reshape(list(NextStep), (len(sequence1)-1 , len(sequence3)-1 ))
    Backward_Directions = nm.vstack([["*"] * Backward_Directions.shape[1], Backward_Directions])
    Backward_Directions = nm.column_stack([["*"] * Backward_Directions.shape[0], Backward_Directions])

    while (True):
        if (Backward_Directions[k][n] == "D"):
            FirstAlign.append(sequence1[k])
            SecondAlign.append(sequence2[k])
            ThirddAlign.append(sequence3[n])
            FourthAlign.append(sequence4[n])
            k -= 1
            n -= 1
        elif (Backward_Directions[k][n] == "U"):
            if (sequence1[k] == ' ' and sequence2[k]==' '):
                FirstAlign.append("-")
                SecondAlign.append("-")
                ThirddAlign.append(sequence3[n])
                FourthAlign.append(sequence4[n])
                n -= 1
            else:
                FirstAlign.append(sequence1[k])
                SecondAlign.append(sequence2[k])
                ThirddAlign.append("-")
                FourthAlign.append("-")
                k -= 1
        else:
            if (sequence3[n] == ' ' and sequence4[n]==' ' ):
                ThirddAlign.append("-")
                FourthAlign.append("-")
                FirstAlign.append(sequence1[k])
                SecondAlign.append(sequence2[k])
                k -= 1
            else:
                ThirddAlign.append(sequence3[n])
                FourthAlign.append(sequence4[n])
                FirstAlign.append("-")
                SecondAlign.append("-")
                n -= 1

        if (k < 0):
            break
        if (n < 0):
            break
    FirstAlign.reverse()
    SecondAlign.reverse()
    ThirddAlign.reverse()
    FourthAlign.reverse()

    print("The Best Alignment Sequences are :")
    print("Alignment Sequence 1 : " , FirstAlign[1:])
    print("Alignment Sequence 2 : " , SecondAlign[1:])
    print("Alignment Sequence 3 : ", ThirddAlign[1:])
    print("Alignment Sequence 4 : ", FourthAlign[1:])
    print("\nThe Scoring Matrix : \n")
    print(Matrix)
    print("The Score Of The Cell We begin From is  : ",  Matrix[len(sequence1) - 1][len(sequence3) - 1], "\n")


       # MAIN FUNCTION
    #-----------------------------------------------------------------------------------------------
print("Please Enter 0 For DNA ALIGNMENT OR 1 for Protein Alignment: ")
Choice = input("Enter Your Choice : ")
if (Choice== '0'):
    print("Enter the four Sequences you want to Align , Please : ")
    Seq1 = input("Enter 1st Sequence : ")
    Seq2 = input("Enter 2nd Sequence : ")
    Seq3 = input("Enter 3nd Sequence : ")
    Seq4 = input("Enter 4nd Sequence : ")
    GAP = input("Enter Gap:")
    print("MSA Of DNA is : \n")
    MSA_DNA(Seq1,Seq2,Seq3,Seq4,GAP)
elif(Choice == '1'):
    print("Enter the four Sequences you want to Align , Please : ")
    Seq1 = input("Enter 1st Sequence : ")
    Seq2 = input("Enter 2nd Sequence : ")
    Seq3 = input("Enter 3nd Sequence : ")
    Seq4 = input("Enter 4nd Sequence : ")
    GAP =input("Enter Gap :")
    print("MSA Of Protein is : \n")
    MSA_protein(Seq1,Seq2,Seq3,Seq4,GAP)
import numpy as nm

from Bio.SubsMat import MatrixInfo
blosum= MatrixInfo.blosum62

def Global_DNA_Align_Algorithm(FirstSeq, SecondSeq):
    # Variables
    Sequence1 = " " + FirstSeq.upper()
    Sequence2 = " " + SecondSeq.upper()
    MScore = 1
    MMScore = -2
    GScore = -1
    FirstAlign = []
    SecondAlign = []
    NextStep = ""
    k = len(Sequence1) - 1
    n = len(Sequence2) - 1
    # Matrix Figuration
    Matrix = nm.zeros((len(Sequence1), len(Sequence2)))
    for i in range(len(Sequence1)):
        Matrix[i][0] = i * GScore
    for j in range(len(Sequence2)):
        Matrix[0][j] = j * GScore

    # Matrix Calculation
    for i in range(1, len(Sequence1)):
        for j in range(1, len(Sequence2)):
            if (Sequence1[i] == Sequence2[j]):
                D = MScore + Matrix[i - 1][j - 1]
            else:
                D = MMScore + Matrix[i - 1][j - 1]
            L = GScore + Matrix[i][j - 1]
            U = GScore + Matrix[i - 1][j]
            MaxScoreIndex = nm.argmax([D, L, U])
            if (MaxScoreIndex == 0):
                NextStep += "D"
            elif (MaxScoreIndex == 1):
                NextStep += "L"
            else:
                NextStep += "U"
            Matrix[i][j] = nm.max([D, L, U])

    # Matrix Traceback
    Backward_Directions = nm.reshape(list(NextStep), (len(Sequence1) - 1, len(Sequence2) - 1))
    Backward_Directions = nm.vstack([["*"] * Backward_Directions.shape[1], Backward_Directions])
    Backward_Directions = nm.column_stack([["*"] * Backward_Directions.shape[0], Backward_Directions])

    while (True):
        if (Backward_Directions[k][n] == "D"):
            FirstAlign.append(Sequence1[k])
            SecondAlign.append(Sequence2[n])
            k -= 1
            n -= 1
        elif (Backward_Directions[k][n] == "U"):
            if (Sequence1[k] == ' '):
                FirstAlign.append("-")
                SecondAlign.append(Sequence2[n])
                n -= 1
            else:
             FirstAlign.append(Sequence1[k])
             SecondAlign.append("-")
             k -= 1
        else:
            if (Sequence2[n] == ' '):
                SecondAlign.append("-")
                FirstAlign.append(Sequence1[k])
                k -= 1
            else:
             SecondAlign.append(Sequence2[n])
             FirstAlign.append("-")
             n -= 1


        if (k < 0):
            break
        if (n < 0):
            break

    FirstAlign.reverse()
    SecondAlign.reverse()
    print("The Best Alignment Sequences are :")
    print("Alignment Sequence 1 : " , FirstAlign[1:])
    print("Alignment Sequence 2 : " , SecondAlign[1:])
    print("\nThe Scoring Matrix : \n")
    print(Matrix)
    print("The Score Of The Cell We begin From is  : ",  Matrix[len(Sequence1) - 1][len(Sequence2) - 1], "\n")






def Global_Protein_Align_Algorithm(FirstSeq, SecondSeq):
        # Variables
        Sequence1 = " " + FirstSeq.upper()
        Sequence2 = " " + SecondSeq.upper()
        GScore = -1
        FirstAlign = []
        SecondAlign = []
        NextStep = ""
        k = len(Sequence1) - 1
        n = len(Sequence2) - 1
        # Matrix Figuration
        Matrix = nm.zeros((len(Sequence1), len(Sequence2)))
        for i in range(len(Sequence1)):
            Matrix[i][0] = i * GScore
        for j in range(len(Sequence2)):
            Matrix[0][j] = j * GScore
        # Matrix Calculation
        for i in range(1, len(Sequence1)):
            for j in range(1, len(Sequence2)):
                pair = (Sequence1[i], Sequence2[j])
                if (Sequence1[i] == Sequence2[j]):
                    if (pair not in blosum):
                        D = blosum[tuple(reversed(pair))] + Matrix[i - 1][j - 1]
                    else:
                        D = blosum[Sequence1[i], Sequence2[j]] + Matrix[i - 1][j - 1]
                else:
                    if (pair not in blosum):
                        D = blosum[tuple(reversed(pair))] + Matrix[i - 1][j - 1]
                    else:
                        D = blosum[Sequence1[i], Sequence2[j]] + Matrix[i - 1][j - 1]
                L = GScore + Matrix[i][j - 1]
                U = GScore + Matrix[i - 1][j]
                MaxScoreIndex = nm.argmax([D, L, U])
                if (MaxScoreIndex == 0):
                    NextStep += "D"
                elif (MaxScoreIndex == 1):
                    NextStep += "L"
                else:
                    NextStep += "U"
                Matrix[i][j] = nm.max([D, L, U])

        # Matrix Traceback
        Backward_Directions = nm.reshape(list(NextStep), (len(Sequence1) - 1, len(Sequence2) - 1))
        Backward_Directions = nm.vstack([["*"] * Backward_Directions.shape[1], Backward_Directions])
        Backward_Directions = nm.column_stack([["*"] * Backward_Directions.shape[0], Backward_Directions])

        while (True):
            if (Backward_Directions[k][n] == "D"):
                FirstAlign.append(Sequence1[k])
                SecondAlign.append(Sequence2[n])
                k -= 1
                n -= 1
            elif (Backward_Directions[k][n] == "U"):
                if (Sequence1[k] == ' '):
                    FirstAlign.append("-")
                    SecondAlign.append(Sequence2[n])
                    n -= 1
                else:
                    FirstAlign.append(Sequence1[k])
                    SecondAlign.append("-")
                    k -= 1
            else:
                if (Sequence2[n] == ' '):
                    SecondAlign.append("-")
                    FirstAlign.append(Sequence1[k])
                    k -= 1
                else:
                    SecondAlign.append(Sequence2[n])
                    FirstAlign.append("-")
                    n -= 1

            if (k < 0):
                break
            if (n < 0):
                break

        FirstAlign.reverse()
        SecondAlign.reverse()
        print("The Best Alignment Sequences are :")
        print("Alignment Sequence 1 : ", FirstAlign[1:])
        print("Alignment Sequence 2 : ", SecondAlign[1:])
        print("\nThe Scoring Matrix : \n")
        print(Matrix)
        print("The Score Of The Cell We begin From is  : ", Matrix[len(Sequence1) - 1][len(Sequence2) - 1], "\n")


        # Main
print("Please Enter 0 For DNA ALIGNMENT OR 1 for Protein Alignment: ")
Choise = input("Enter Your Choise : ")
if (Choise=='0'):
    print("Enter the Two Sequences you want to Align , Please : ")
    Seq1 = input("Enter 1st Sequence : ")
    Seq2 = input("Enter 2nd Sequence : ")
    print("Global Alignment is : \n")
    Global_DNA_Align_Algorithm(Seq1, Seq2)
elif(Choise=='1'):
    print("Enter the Two Sequences you want to Align , Please : ")
    Seq1 = input("Enter 1st Sequence : ")
    Seq2 = input("Enter 2nd Sequence : ")
    print("Global Alignment is : \n")
    Global_Protein_Align_Algorithm(Seq1, Seq2)


import numpy as nm
from Bio.SubsMat import MatrixInfo
blosum= MatrixInfo.blosum62

def Local_DNA_Align_Algorithm(FirstSeq, SecondSeq):
    # Variables
    Seq1 = " " +FirstSeq.upper()
    Seq2 = " " +SecondSeq.upper()
    MScore = 1
    MMScore = -2
    GScore = -1
    FirstAlign = []
    SecondAlign = []
    NextStep = ""
    # Matrix Figuration
    Matrix = nm.zeros((len(Seq1), len(Seq2)))
    for i in range(1, len(Seq1)):
        for j in range(1, len(Seq2)):
            if (Seq1[i] == Seq2[j]):
                D = MScore + Matrix[i - 1][j - 1]
            else:
                D = MMScore + Matrix[i - 1][j - 1]
            L = GScore + Matrix[i][j - 1]
            U = GScore + Matrix[i - 1][j]
            MaxScoreIndex = nm.argmax([D, L, U])  # Index of First Max Score of the Current Step
            if (MaxScoreIndex == 0):
                NextStep += "D"
            elif (MaxScoreIndex == 1):
                NextStep += "L"
            else:
                NextStep += "U"

            Matrix[i][j] = nm.max([D, L, U])  # Max Score of the matrix

            if (Matrix[i][j] < 0):
                Matrix[i][j] = 0

    # Matrix Traceback

    Backward_Directions = nm.reshape(list(NextStep), (len(Seq1) - 1, len(Seq2) - 1 ))
    Backward_Directions = nm.vstack([["0"] * Backward_Directions.shape[1], Backward_Directions])
    Backward_Directions = nm.column_stack([["0"] * Backward_Directions.shape[0], Backward_Directions])

    #Getting the index of the Highest Score

    MaxValue = nm.max(Matrix)
    for k in range(1, len(Seq1)):
        for n in range(1, len(Seq2)):
            if (Matrix[k][n] == MaxValue):
                y = k
                z = n

    #Getting the Alingmnet Sequences

    while (True):
        if (Backward_Directions[y][z] == "D"):
            FirstAlign.append(Seq1[y])
            SecondAlign.append(Seq2[z])
            y -= 1
            z -= 1
        elif (Backward_Directions[y][z] == "U"):
            FirstAlign.append(Seq1[y])
            SecondAlign.append("-")
            y -= 1
        else:
            SecondAlign.append(Seq2[z])
            FirstAlign.append("-")
            z -=1

        if (y < 0):
            break
        if (z < 0):
            break
        if (Matrix[y][z] == 0):
            break

    FirstAlign.reverse()
    SecondAlign.reverse()


    print("The Best Matched Sequences are :")
    print("First Sequence : " ,FirstAlign)
    print("Second Sequence : " ,SecondAlign)
    print("\nThe Scoring Matrix : \n")
    print(Matrix)
    print("Optimal Alignment Score Where Traceback From : ", MaxValue, "\n")

def Local_Protein_Align_Algorithm(FirstSeq, SecondSeq):

        # Variables
        Seq1 = " " + FirstSeq.upper()
        Seq2 = " " + SecondSeq.upper()
        GScore = -1
        FirstAlign = []
        SecondAlign = []
        NextStep = ""
        # Matrix Figuration
        Matrix = nm.zeros((len(Seq1), len(Seq2)))
        for i in range(1, len(Seq1)):
            for j in range(1, len(Seq2)):
                pair = (Seq1[i], Seq2[j])
                if (Seq1[i] == Seq2[j]):
                    if(pair not in blosum):
                        D = blosum[tuple(reversed(pair))]+Matrix[i-1][j-1]
                    else:
                        D = blosum[Seq1[i],Seq2[j]]+Matrix[i-1][j-1]
                else:
                    if (pair not in blosum):
                        D = blosum[tuple(reversed(pair))] + Matrix[i - 1][j - 1]
                    else:
                        D = blosum[Seq1[i], Seq2[j]] + Matrix[i - 1][j - 1]
                L = GScore + Matrix[i][j - 1]
                U = GScore + Matrix[i - 1][j]
                MaxScoreIndex = nm.argmax([D, L, U])  # Index of First Max Score of the Current Step
                if (MaxScoreIndex == 0):
                    NextStep += "D"
                elif (MaxScoreIndex == 1):
                    NextStep += "L"
                else:
                    NextStep += "U"

                Matrix[i][j] = nm.max([D, L, U])  # Max Score of the matrix

                if (Matrix[i][j] < 0):
                    Matrix[i][j] = 0

        # Matrix Traceback

        Backward_Directions = nm.reshape(list(NextStep), (len(Seq1) - 1, len(Seq2) - 1))
        Backward_Directions = nm.vstack([["0"] * Backward_Directions.shape[1], Backward_Directions])
        Backward_Directions = nm.column_stack([["0"] * Backward_Directions.shape[0], Backward_Directions])

        # Getting the index of the Highest Score

        MaxValue = nm.max(Matrix)
        for k in range(1, len(Seq1)):
            for n in range(1, len(Seq2)):
                if (Matrix[k][n] == MaxValue):
                    y = k
                    z = n

        # Getting the Alingmnet Sequences

        while (True):
            if (Backward_Directions[y][z] == "D"):
                FirstAlign.append(Seq1[y])
                SecondAlign.append(Seq2[z])
                y -= 1
                z -= 1

            elif (Backward_Directions[y][z] == "U"):
                FirstAlign.append(Seq1[y])
                SecondAlign.append("-")
                y -= 1
            else:
                SecondAlign.append(Seq2[z])
                FirstAlign.append("-")
                z -= 1

            if (y < 0):
                break
            if (z < 0):
                break
            if (Matrix[y][z] == 0):
                break

        FirstAlign.reverse()
        SecondAlign.reverse()

        print("The Best Matched Sequences are :")
        print("First Sequence : ", FirstAlign)
        print("Second Sequence : ", SecondAlign)
        print("\nThe Scoring Matrix : \n")
        print(Matrix)
        print("Optimal Alignment Score Where Traceback From : ", MaxValue, "\n")

    # Main


print("Enter 0 For DNA Alignment OR 1 For Protein Alignment : ")
Choice = input("Enter Your Choice : " )
print("Enter the Two Sequences you want to Align , Please : ")
Seq1 = input("Enter 1st Sequence : ")
Seq2 = input("Enter 2nd Sequence : ")
if (Choice=='0'):
    print("Local Alignment For DNA is : \n")
    Local_DNA_Align_Algorithm(Seq1,Seq2)
elif(Choice=='1'):
    print("Local Alignment For Protein is  : \n")
    Local_Protein_Align_Algorithm(Seq1,Seq2)

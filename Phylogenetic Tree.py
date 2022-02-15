def Alphabetical_Order(Start, End):
    ClusterPairs = []
    for i in range(ord(Start), ord(End) + 1):  # Ord Fun : Take character argument and returns the Unicode Code of it
        ClusterPairs.append(
            chr(i))  # Chr Fun : Take Integer argument (Unicode Value) and returns Character representing its value
    return ClusterPairs


def MaTMiniCell(Matrix):
    # infinite float -> useful for finding lowest values for something like in matrix & handle -ve , +ve float numbers
    MinCell = float("inf")
    a = -1
    b = -1
    for i in range(len(Matrix)):
        for j in range(len(Matrix[i])):
            if Matrix[i][j] < MinCell:
                MinCell = Matrix[i][j]
                a = i
                b = j

    return a, b  # Return the indexes of the First Min Cell in the Matrix


def PairingOfLeaves(ClusterPairs, F_Index, S_Index):
    if S_Index < F_Index:
        F_Index, S_Index = S_Index, F_Index  # The Way I Swap the Indexes In Python

    # Combining the Pairs In 1st Index of the Paired Leaves
    ClusterPairs[F_Index] = "(" + ClusterPairs[F_Index] + " , " + ClusterPairs[S_Index] + ")"

    del ClusterPairs[S_Index]  # Remove the Immediately redundant pair Part in the 2nd Index

    print(ClusterPairs)


def RebuildingZMatrix(matrix, F_Index, S_Index):
    if S_Index < F_Index:
        F_Index, S_Index = S_Index, F_Index  # Swapping the Indexes
    row = []
    print("Distance Matrix : " ,matrix, "\n")
    # Rebuild the Row indexes
    for i in range(0, F_Index):

        # Calculate Avg. of two cell that they in the same column and after finishing this we del the redundant row at last step
        row.append((matrix[F_Index][i] + matrix[S_Index][i]) / 2)

    # Put the Paired Leaves in the 1st Index of row  ex. if f.i=0 , s.i=1  --> index of combined leaves is = 0
    matrix[F_Index] = row

    # Hint : Second row now contains the values of the row below and so on  Where 2nd_Row_indexes < S_index for Column
    for i in range(F_Index + 1, S_Index):

        matrix[i][F_Index] = (matrix[i][F_Index] + matrix[S_Index][i]) / 2

    # Calculate the remaining rows' Values given the row position of the min cell
    for i in range(S_Index + 1, len(matrix)):

        matrix[i][F_Index] = (matrix[i][F_Index] + matrix[i][S_Index]) / 2

        del matrix[i][S_Index]  # Remove Immediately the redundant 2nd Index Column


    del matrix[S_Index]  # Remove Immediately the redundant 2nd Index row



def UPGMA(matrix, ClusterPairs):
    while ( len(ClusterPairs) > 1 ):
        i, j = MaTMiniCell(matrix)
        PairingOfLeaves(ClusterPairs, i, j)  # Update Cluster Pairs
        RebuildingZMatrix(matrix, i, j)  # Combine and rebuild Matrix Values
    print("\nFinal Clusters : ")
    return ClusterPairs[0]  # Return the final Cluster Pairs

    # Main


Matrix_Characters = Alphabetical_Order("A", "E")
MATRIX = [
    [],  # A
    [20],  # B
    [60, 50],  # C
    [100, 90, 40],  # D
    [90, 80, 50, 30],  # E
]
print(UPGMA(MATRIX, Matrix_Characters))

''' []   j=0  j=1  j=2  j=3  j=4
     -    A    B    C    D    E
i=0  A    0    

i=1  B    20   0

i=2  C    60   50   0

i=3  D    100  90  40    0

i=4  E    90   80   50   30   0

Length of Matrix is = 5


'''

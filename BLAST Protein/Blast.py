################################## Main ##################################
query_seq =input("Enter Protein Query Sequence : ")                            # "PQGEFG"
len_of_word = input("Enter The Length Of Word : ")                         #4
word_threshold = input("Enter Word Threshold  : ")                      #13
HSP_threshold = input("Enter HSP Threshold : ")                       #10

def RemoveRepeatedSeq(query_seq):
    query_seq = query_seq.upper()
    query_list = list(query_seq)
    firstIndex = 0
    c = 1
    for i in range(1, len(query_list)):
        if (query_list[i] == query_list[i - 1]):
            c = c + 1
        if (query_list[i] != query_list[i - 1] or i == len(query_list) - 1):
            if (c >= 4):
                if (i == len(query_list) - 1 and query_list[i] == query_list[i - 1]):
                    i = i + 1
                for j in range(firstIndex, i):
                    query_list[j] = "X"
            c = 1
            firstIndex = i
    b = ''.join(query_list)  # make a list of strings in a single string
    # final = []
    # for i in range(0, len(query_list)):
    # if (query_list[i] != "X"):
    #  final.append(query_list[i])
    #   b= ''.join(final) # make a ilst of strings in a single string
    return b


# query_seq = "TCAPNGGGGARGKLVGGGGMMMENNNN"   #input
# query_sequence = 'AGTNNNNNTAGAC'
queryRes = RemoveRepeatedSeq(query_seq)
print(queryRes)


#len_of_word = 3   #K    # input
def MaKeWords(queryRes):
	words = []
	N = len(queryRes)
	for i in range (0,N-int(len_of_word)+1):   # num of words = n - K + 1
		if (queryRes[i:i+int(len_of_word)] not in words):
			words.append(queryRes[i:i+int(len_of_word)])
	return words
words=MaKeWords(queryRes)
print(words)

AminoAcids = ["A", "R", "N", "D", "C", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "E"]


def getPossibleCombination(words):
    listoflists = []
    for word in words:

        combination = []
        combination.insert(0, word)  # maKe the orginal word the first in the list
        for i in range(0, int(len_of_word)):
            for aminoacid in AminoAcids:

                if (word[i] == aminoacid and word[i] not in combination):

                    combination.append(word)
                else:
                    char = word
                    char = char[:i] + aminoacid + char[i + 1:]  # if i =0 char[:i]=0
                    combination.append(char)
        listoflists.append((list(combination)))
    return listoflists

neiborhood_words = getPossibleCombination(words)
print(neiborhood_words)



# scoring matrix from ncbi
blosum62 = {
    '*':{'*':1,'A':-4,'C':-4,'B':-4,'E':-4,'D':-4,'G':-4,'F':-4,'I':-4,'H':-4,'K':-4,'M':-4,'L':-4,'N':-4,'Q':-4,'P':-4,'S':-4,'R':-4,'T':-4,'W':-4,'V':-4,'Y':-4,'X':-4,'Z':-4},
    'A':{'*':-4,'A':4,'C':0,'B':-2,'E':-1,'D':-2,'G':0,'F':-2,'I':-1,'H':-2,'K':-1,'M':-1,'L':-1,'N':-2,'Q':-1,'P':-1,'S':1,'R':-1,'T':0,'W':-3,'V':0,'Y':-2,'X':-1,'Z':-1},
    'C':{'*':-4,'A':0,'C':9,'B':-3,'E':-4,'D':-3,'G':-3,'F':-2,'I':-1,'H':-3,'K':-3,'M':-1,'L':-1,'N':-3,'Q':-3,'P':-3,'S':-1,'R':-3,'T':-1,'W':-2,'V':-1,'Y':-2,'X':-1,'Z':-3},
    'B':{'*':-4,'A':-2,'C':-3,'B':4,'E':1,'D':4,'G':-1,'F':-3,'I':-3,'H':0,'K':0,'M':-3,'L':-4,'N':3,'Q':0,'P':-2,'S':0,'R':-1,'T':-1,'W':-4,'V':-3,'Y':-3,'X':-1,'Z':1},
    'E':{'*':-4,'A':-1,'C':-4,'B':1,'E':5,'D':2,'G':-2,'F':-3,'I':-3,'H':0,'K':1,'M':-2,'L':-3,'N':0,'Q':2,'P':-1,'S':0,'R':0,'T':-1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':4},
    'D':{'*':-4,'A':-2,'C':-3,'B':4,'E':2,'D':6,'G':-1,'F':-3,'I':-3,'H':-1,'K':-1,'M':-3,'L':-4,'N':1,'Q':0,'P':-1,'S':0,'R':-2,'T':-1,'W':-4,'V':-3,'Y':-3,'X':-1,'Z':1},
    'G':{'*':-4,'A':0,'C':-3,'B':-1,'E':-2,'D':-1,'G':6,'F':-3,'I':-4,'H':-2,'K':-2,'M':-3,'L':-4,'N':0,'Q':-2,'P':-2,'S':0,'R':-2,'T':-2,'W':-2,'V':-3,'Y':-3,'X':-1,'Z':-2},
    'F':{'*':-4,'A':-2,'C':-2,'B':-3,'E':-3,'D':-3,'G':-3,'F':6,'I':0,'H':-1,'K':-3,'M':0,'L':0,'N':-3,'Q':-3,'P':-4,'S':-2,'R':-3,'T':-2,'W':1,'V':-1,'Y':3,'X':-1,'Z':-3},
    'I':{'*':-4,'A':-1,'C':-1,'B':-3,'E':-3,'D':-3,'G':-4,'F':0,'I':4,'H':-3,'K':-3,'M':1,'L':2,'N':-3,'Q':-3,'P':-3,'S':-2,'R':-3,'T':-1,'W':-3,'V':3,'Y':-1,'X':-1,'Z':-3},
    'H':{'*':-4,'A':-2,'C':-3,'B':0,'E':0,'D':-1,'G':-2,'F':-1,'I':-3,'H':8,'K':-1,'M':-2,'L':-3,'N':1,'Q':0,'P':-2,'S':-1,'R':0,'T':-2,'W':-2,'V':-3,'Y':2,'X':-1,'Z':0},
    'K':{'*':-4,'A':-1,'C':-3,'B':0,'E':1,'D':-1,'G':-2,'F':-3,'I':-3,'H':-1,'K':5,'M':-1,'L':-2,'N':0,'Q':1,'P':-1,'S':0,'R':2,'T':-1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':1},
    'M':{'*':-4,'A':-1,'C':-1,'B':-3,'E':-2,'D':-3,'G':-3,'F':0,'I':1,'H':-2,'K':-1,'M':5,'L':2,'N':-2,'Q':0,'P':-2,'S':-1,'R':-1,'T':-1,'W':-1,'V':1,'Y':-1,'X':-1,'Z':-1},
    'L':{'*':-4,'A':-1,'C':-1,'B':-4,'E':-3,'D':-4,'G':-4,'F':0,'I':2,'H':-3,'K':-2,'M':2,'L':4,'N':-3,'Q':-2,'P':-3,'S':-2,'R':-2,'T':-1,'W':-2,'V':1,'Y':-1,'X':-1,'Z':-3},
    'N':{'*':-4,'A':-2,'C':-3,'B':3,'E':0,'D':1,'G':0,'F':-3,'I':-3,'H':1,'K':0,'M':-2,'L':-3,'N':6,'Q':0,'P':-2,'S':1,'R':0,'T':0,'W':-4,'V':-3,'Y':-2,'X':-1,'Z':0},
    'Q':{'*':-4,'A':-1,'C':-3,'B':0,'E':2,'D':0,'G':-2,'F':-3,'I':-3,'H':0,'K':1,'M':0,'L':-2,'N':0,'Q':5,'P':-1,'S':0,'R':1,'T':-1,'W':-2,'V':-2,'Y':-1,'X':-1,'Z':3},
    'P':{'*':-4,'A':-1,'C':-3,'B':-2,'E':-1,'D':-1,'G':-2,'F':-4,'I':-3,'H':-2,'K':-1,'M':-2,'L':-3,'N':-2,'Q':-1,'P':7,'S':-1,'R':-2,'T':-1,'W':-4,'V':-2,'Y':-3,'X':-1,'Z':-1},
    'S':{'*':-4,'A':1,'C':-1,'B':0,'E':0,'D':0,'G':0,'F':-2,'I':-2,'H':-1,'K':0,'M':-1,'L':-2,'N':1,'Q':0,'P':-1,'S':4,'R':-1,'T':1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':0},
    'R':{'*':-4,'A':-1,'C':-3,'B':-1,'E':0,'D':-2,'G':-2,'F':-3,'I':-3,'H':0,'K':2,'M':-1,'L':-2,'N':0,'Q':1,'P':-2,'S':-1,'R':5,'T':-1,'W':-3,'V':-3,'Y':-2,'X':-1,'Z':0},
    'T':{'*':-4,'A':0,'C':-1,'B':-1,'E':-1,'D':-1,'G':-2,'F':-2,'I':-1,'H':-2,'K':-1,'M':-1,'L':-1,'N':0,'Q':-1,'P':-1,'S':1,'R':-1,'T':5,'W':-2,'V':0,'Y':-2,'X':-1,'Z':-1},
    'W':{'*':-4,'A':-3,'C':-2,'B':-4,'E':-3,'D':-4,'G':-2,'F':1,'I':-3,'H':-2,'K':-3,'M':-1,'L':-2,'N':-4,'Q':-2,'P':-4,'S':-3,'R':-3,'T':-2,'W':11,'V':-3,'Y':2,'X':-1,'Z':-3},
    'V':{'*':-4,'A':0,'C':-1,'B':-3,'E':-2,'D':-3,'G':-3,'F':-1,'I':3,'H':-3,'K':-2,'M':1,'L':1,'N':-3,'Q':-2,'P':-2,'S':-2,'R':-3,'T':0,'W':-3,'V':4,'Y':-1,'X':-1,'Z':-2},
    'Y':{'*':-4,'A':-2,'C':-2,'B':-3,'E':-2,'D':-3,'G':-3,'F':3,'I':-1,'H':2,'K':-2,'M':-1,'L':-1,'N':-2,'Q':-1,'P':-3,'S':-2,'R':-2,'T':-2,'W':2,'V':-1,'Y':7,'X':-1,'Z':-2},
    'X':{'*':-4,'A':-1,'C':-1,'B':-1,'E':-1,'D':-1,'G':-1,'F':-1,'I':-1,'H':-1,'K':-1,'M':-1,'L':-1,'N':-1,'Q':-1,'P':-1,'S':-1,'R':-1,'T':-1,'W':-1,'V':-1,'Y':-1,'X':-1,'Z':-1},
    'Z':{'*':-4,'A':-1,'C':-3,'B':1,'E':4,'D':1,'G':-2,'F':-3,'I':-3,'H':0,'K':1,'M':-1,'L':-3,'N':0,'Q':3,'P':-1,'S':0,'R':0,'T':-1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':4}}
final_words = []
scores = []
#all_scores = []
#word_threshold = 13 # input
def ComputeScores(neiborhood_words):

    for list in range (0,len(neiborhood_words)):
        for i in range(0,len(neiborhood_words[list])):
            total = 0
            for j in range(0,int(len_of_word)):
                    total += blosum62[neiborhood_words[list][0][j]][ neiborhood_words[list][i][j]]
               # except KeyError:
                #    total += blosum62[(neiborhood_words[list][i][j], neiborhood_words[list][0][j])]
                    if (total >=int(int(word_threshold))):      # for int(word_threshold)  , thershold =10
                      scores.append(total)  # for thershold
                      final_words.append(neiborhood_words[list][i])  # for int(word_threshold)
         #   all_scores.append(total)    # h7tgha fel extend
            print((neiborhood_words[list][i] , total))   # print all scores
            dic = dict(zip(final_words, scores)) #for thershold
    return dic
finalss =ComputeScores(neiborhood_words)
print ("Final Seeds :  " + str(finalss))
#print(ComputeScores(neiborhood_words))


#open database from file
file = open("database_Lab.txt", "r")
#file = open("Protein_database.txt", "r")
list = [(line.strip()).split() for line in file]
file.close()
targetSeq = []
for i in range(0, len(list)):
    targetSeq.append(list[i][0])

'''targetSeq = [
 "KMMKSFFLVVTILALTLPFLGAQEQNQEQPIRCEKDERFFSDKIAKYIPIQYVLSRYPSYGLNYYQQKPVALINNQFLPYPYYAKPAAVRSPAQILQWQVLSNTVPAKSCQAQPTTMARHPHPHLSFMAIPPKKNQDKTEIPTINTIASGEPTSTPTTEAVESTVATLEDSPEVIESPPEINTVQVTSTAV",
 "MKLFWLLFTIGFCWAQYSSNTQQGRTSIVHLFEWRWVDIALECERYLAPKGFGGVQVSPPNENVAIHNPFRPWWERYQPVSYKLCTRSGNEDEFRNMVTRCNNVGVRIYVDAVINHMCGNAVSAGTSSTCGSYFNPGSRDFPAVPYSGWDFNDGKCKTGSGDIENYNDATQVRDCRLSGLLDLALGKDYVRSKIAEYMNHLIDIGVAGFRIDASKHMWPGDIKAILDKLHNLNSNWFPEGSKPFIYQEVIDLGGEPIKSSDYFGNGRVTEFKYGAKLGTVIRKWNGEKMSYLKNWGEGWGFMPSDRALVFVDNHDNQRGHGAGGASILTFWDARLYKMAVGFMLAHPYGFTRVMSSYRWPRYFENGKDVNDWVGPPNDNGVTKEVTINPDTTCGNDWVCEHRWRQIRNMVNFRNVVDGQPFTNWYDNGSNQVAFGRGNRGFIVFNNDDWTFSLTLQTGLPAGTYCDVISGDKINGNCTGIKIYVSDDGKAHFSISNSAED",
 "MKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNELSKDIGSESTEDQAMEDIKQMEAESISSSEEIVPNSVEQKHIQKEDVPSERYLGYLEQLLRLKKYKVPQLEIVPNSAEERLHSMKEGIHAQQKEPMIGVNQELAYFYPELFRQFYQLDAYPSGAWYYVPLGTQYTDAPSFSDIPNPIGSENSEKTTMPLW" ,
 "MKFFIFTCLLAVALAKNTMEHVSSSEESIISQETYKQEKNMAINPSKENLCSTFCKEVVRNANEEEYSIGSSSEESAEVATEEVKITVDDKHYQKALNEINQFYQKFPQYLQYLYQGPIVLNPWDQVKRNAVPITPTLNREQLSTSEENSKKTVDMESTEVFTKKTKLTEEEKNRLNFLKKISQRYQKFALPQYLKTVYQHQKAMKPWIQPKTKVIPYVRYL" ,
 ]'''

target_info = []
firstindex = []
lastindex = []
seedHit = []
numseq = []

def checKHit(targetSeq, final_words):
    for word in final_words:
        for seq in targetSeq:
            for i in range(0, len(seq)):
                if (word == seq[i:i + int(len_of_word)] and word not in seedHit):
                    seedHit.append(word)
                    numseq.append(seq)
                    dfirst = i
                    dlast = i + (int(len_of_word) - 1)
                    firstindex.append(i)
                    lastindex.append(i + (int(len_of_word) - 1))
                    ID = targetSeq.index(seq)
                    target_info.append([word, seq, dfirst, dlast, ID])
            dic = dict(zip(seedHit, numseq))  # hit and which seq
            Kic = dict(zip(firstindex, lastindex))  # first and last index for each hit

    return dic, Kic


print(checKHit(targetSeq, final_words))
if (len(seedHit) == 0):
    print("The Query Was Not Found in Database !")

query_info = []
query_info1 = []

for seed in seedHit:
    for list in range(0, len(neiborhood_words)):
        for i in range(0, len(neiborhood_words[list])):
            if (seed == neiborhood_words[list][i]):
                original_word = neiborhood_words[list][0]

                for j in range(0, len(queryRes)):
                    if (queryRes[j:j + int(len_of_word)] == original_word):
                        original_firstindex = j
                        original_lastindex = j + (int(len_of_word) - 1)
                        if (seed not in query_info1):  # and seed not in query_info[i][i]):
                            query_info1.append([seed, original_word, queryRes, original_firstindex, original_lastindex])
for each in query_info1:
    if each not in query_info:
        query_info.append(each)

# print(t_score)
print(query_info)
print(target_info)
# print (blosum62['P']['S'])
# print(target_info[0])


def extend(seedHit, target_info, query_info):
    for seed in range(0, len(seedHit)):
        q_score = 0       # seed score
        RChar_score = 0   #right char score
        Char_score = 0    #left char score

        # HSP_threshold=13  # input
        new_queryfirstindex = query_info[seed][3] - 1
        new_querylastindex = query_info[seed][4] + 1
        new_targetfirstindex = target_info[seed][2] - 1
        new_targetlastindex = target_info[seed][3] + 1
        Result = [query_info[seed][1]]
        T_Result = [query_info[seed][0]]

        result = ""
        T_result = ""

        for i in range(0, int(len_of_word)):
            q_score += blosum62[query_info[seed][0][i]][query_info[seed][1][i]]
            Final_score = q_score

        HSP = Final_score
        for j in range(0, len(query_info[seed][2])):
                if (target_info[seed][2] != 0 and query_info[seed][3] != 0 and new_queryfirstindex > -1 and new_targetfirstindex > -1):  # extend from left  query
                    q_leftchar = query_info[seed][2][new_queryfirstindex]
                    t_leftchar = target_info[seed][1][new_targetfirstindex]

                    Char_score = blosum62[q_leftchar][t_leftchar]

                    if (HSP - (Final_score + Char_score) <= 3):
                        Final_score = Char_score + Final_score
                        if (Final_score > HSP):
                            HSP = Final_score
                        Result.insert(0, q_leftchar)
                        T_Result.insert(0, t_leftchar)
                        new_queryfirstindex -= 1
                        new_targetfirstindex -= 1

                        print("Score : " + str(Final_score) + " Left HSP : " + str(HSP) + "  HSP - FINAL : " + str(HSP - Final_score))

                if (target_info[seed][3] != len(target_info[seed][1]) and query_info[seed][4] != len( query_info[seed][2]) and new_querylastindex < len(query_info[seed][2]) and new_targetlastindex < len(target_info[seed][1])):  # extend from right  query
                    q_rightchar = query_info[seed][2][new_querylastindex]
                    t_rightchar = target_info[seed][1][new_targetlastindex]
                    RChar_score = blosum62[q_rightchar][t_rightchar]

                    if (HSP - (Final_score + RChar_score) <= 3):
                        new_querylastindex += 1
                        new_targetlastindex += 1

                        Final_score = Final_score + RChar_score
                        print("Score : " + str(Final_score) + " Right HSP : " + str(HSP) + "  HSP - FINAL : " + str(HSP - Final_score))
                        if (Final_score > HSP):
                            HSP = Final_score

                        Result.append(q_rightchar)
                        T_Result.append(t_rightchar)

        if (Final_score >= int(HSP_threshold)):
            print("Seed Score      : " + str(q_score))
            print("Query Sequance  : " + result.join(Result))
            print("Target Sequance : " + T_result.join(T_Result))
            print("Sequance ID     : " + str(target_info[seed][4]))
            print("Final_score     : " + str(Final_score))
            print("************************")


K = extend(seedHit, target_info, query_info)
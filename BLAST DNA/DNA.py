#na2s awl step ale bshel feha low complexity regions + get seeds remove which its score < S + extend and get score
Match ={
    'A':{'A':1, 'T':2, 'G':3, 'C':4},
    'T':{'A':2, 'T':1, 'G':2, 'C':2},
    'G':{'A':3, 'T':2, 'G':1, 'C':3},
    'C':{'A':4, 'T':2, 'G':3, 'C':1}
}
#multiple sequance every sequance in line
def Database(filename):
    file=open(filename)
    DB=[]
    for line in file:
        DB.append(line.rstrip())
    file.close()
    return DB

#one sequance in one line
def Sequance(filename):
    file = open(filename)
    Seq=(file.read())
    file.close()
    return Seq

def GetWords(seq):
    Words=[]
    for i in range(len(seq)):
        if len(seq[i:i+11]) < 11:
            continue
        else:
            Words.append(seq[i:i+11])
    return Words

def Checkwords(Words,Seq):
    if len(Words)==len(Seq)-11+1:
        return True
    else:
        return False

def GetNeighborhood(Words,T):
    list=['A','T','C','G']
    Neighborhood=[]
    seeds = []
    for i in range (len(Words)):
        for j in range (len(Words[i])):
            for k in range(len(list)):
                score = 0
                add=Words[i].replace(Words[i][j],list[k])
                if len(add)==11 and add not in Neighborhood:
                    for x in range (len(add)):
                        score+=Match[add[x]][Words[i][x]]
                    #print("word",Words[i],"add",add,"score",score)
                    Neighborhood.append(add)
                    if score >= T :
                        seeds.append((add,i))
    return Neighborhood,seeds

#hit contain [0]index in query [1]which seq in DB [2]which index in DB seq [3]seed ale 3ml hit
def Gethits(Db,seeds):
    Hitindex=[]
    for i in range(len(seeds)):
        for j in range(len(Db)):
            if seeds[i][0] in Db[j]:
                indexindb=Db[j].find(seeds[i][0])
                Hitindex.append((seeds[i][1],j,indexindb,seeds[i][0]))
    return Hitindex

def Extend(hits,seq,db):
    HSP = []
    #################Forword#################
    for i in range (len(hits)):
        finalSEQ = ""
        finalDb = ""
        count = 1
        oldscore = 0
        newscore = 0
        forword = 0
        backword = 0
        #hmshy forword w backword l7d emta
        maxcount = 0
        if len(seq) < len(db):
            minlen = len(seq)
            forword = minlen - 11 - (hits[i][0])
            backword = len(seq) - 11 - forword
            maxcount = min(forword,backword)
        elif len(db[hits[i][1]]) < len(seq):
            minlen = len(db[hits[i][1]])
            forword = minlen - 11 - (hits[i][2])
            backword = len(db[hits[i][1]]) - 11 - forword
            maxcount = min(forword,backword)


        #old score between seed and hit in db
        sindex = hits[i][0]
        dseq = db[hits[i][1]]
        dindex = hits[i][2]
        query = hits[i][3]
        dbbbb = dseq[dindex:dindex + 11]
        finalSEQ += query
        finalDb += dbbbb
        for k in range(len(hits[i][3])):
            oldscore += Match[query[k]][dbbbb[k]]

        #can't extend in 2 direction
        if maxcount == 0:
            #no extend
            if forword == 0 and backword == 0:
                HSP.append((query,dbbbb,oldscore,"noextend"))
            #extend forword
            elif forword > 0 and backword == 0:
                countf = 0
                newscore += oldscore
                for a in range(forword):
                    # 7arf 7arf l2odam
                    queryword = seq[sindex + 11 + countf]
                    ##
                    dbword = dseq[dindex + 11 + countf]

                    # newscore
                    newscore += Match[queryword][dbword]
                    if newscore < oldscore:
                        HSP.append((finalSEQ, finalDb, oldscore,"f"))
                        continue

                    else:
                        finalSEQ += queryword
                        finalDb += dbword
                        oldscore = newscore

                    countf += 1
                HSP.append((finalSEQ, finalDb, newscore,"f"))
            #extend backword
            elif forword == 0 and backword > 0:
                countb = 0
                newscore += oldscore
                for j in range(backword):
                    # 7arf 7arf lwara
                    queryword = seq[sindex - countb]
                    ##
                    dbword = dseq[dindex - countb]
                    #newscore
                    newscore += Match[queryword][dbword]
                    if newscore < oldscore:
                        finalSEQ += query
                        finalDb += dbbbb
                        HSP.append((finalSEQ, finalDb, newscore, "b"))
                        continue

                    else:
                        finalSEQ += queryword
                        finalDb += dbword
                        oldscore = newscore

                    countb += 1


        #extend in  2 direction
        elif maxcount>0:
            for j in range(maxcount):
                #forword&backword
                sindex = hits[i][0]
                #7tet 7arf mn wara + seed + 7arf mn 2odam
                queryword = seq[sindex - count]
                queryword += query
                queryword += seq[sindex+2+count]
                ##
                dbword = dseq[dindex - count]
                dbword += dbbbb
                dbword += dseq[dindex+2+count]

                #newscore
                for l in range(len(dbword)):
                    newscore += Match[queryword[l]][dbword[l]]
                if newscore < oldscore:
                    HSP.append((finalSEQ,finalDb,oldscore,"f&b"))
                    continue

                else:
                    HSP.append((queryword,dbword,newscore,"f&b"))

                count += 1
    return HSP

def CalcScore(H):
    Score = 0
    for i in range (len(H)):
        Score += H[i][2]
    return Score



Db=Database("DNADatabase.txt")
Seq=Sequance("DNA.txt")
Words=GetWords(Seq)
Check = Checkwords(Words,Seq)
neighborhood,seeds=GetNeighborhood(Words,1)
hits=Gethits(Db,seeds)
Hsps = Extend(hits,Seq,Db)
print(hits)
print(Hsps)
print("Score =",CalcScore(Hsps))
'''print(data)
print(Seq)
print(Words)
print(Check)
print(neighborhood)
print(seeds)
#print(len(neighborhood))
print(hits)'''
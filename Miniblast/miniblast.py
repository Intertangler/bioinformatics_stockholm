#first step is to create a database
#the database is a simple list of sequences (strings)
#we define a function that loads the list from a text file
def read_database(filename):
    f = open(filename)
    db = []
    for line in f:
        db.append(line.rstrip())
    f.close()
    return db

#next function creates a dictionary - a map that will identify all words of size w
#we keep track of the positions where it occurs
#this technique is known as hashing
def build_map(query, w):
    res = {}
    for i in range(len(query) - w + 1):
        subseq = query[i:i+w]
        if subseq in res:
            res[subseq].append(i)
        else:
            res[subseq] = [i]
    return res

#in our simplified version of blast, we will not use a substitution matrix
#we will not consider gaps
#only scores of 1 or 0, match of mismatch will be considered
def get_hits(seq,m,w):
    res = []
    for i in range(len(seq)-w+1):
        subseq = seq[i:i+w]
        if subseq in m:
            l = m[subseq]
            for ind in l:
                res.append((ind,i))
    return res
def extends_hit ( seq, hit, query, w):
    stq, sts = hit[0], hit[1]
    matfw = 0
    k = 0
    bestk = 0
    while 2*matfw >= k and stq+w+k < len(query) and sts+w+k < len(seq):
        if query[stq+w+k] == seq[sts+w+k]:
            matfw+=1
            bestk=k+1
        k+= 1
    size = w + bestk
    k = 0
    matbw = 0
    bestk = 0
    while 2*matbw >= k and stq > k and sts > k:
        if query[stq-k-1] == seq[sts-k-1]:
            matbw += 1
            bestk = k +1
        k+=1
    size += bestk
    return (stq-bestk,sts-bestk,size,w+matfw+matbw)
def hit_best_score(seq,query,m,w):
    hits = get_hits(seq,m,w)
    bestScore = -1.0
    best = ()
    for h in hits:
        ext = extends_hit(seq,h,query,w)
        score = ext[3]
        if score > bestScore or (score == bestScore and ext[2] < best[2]):
            bestScore = score
            best = ext
    return best
def best_alignment(db,query,w):
    m = build_map(query,w)
    bestScore = -1.0
    res = (0,0,0,0,0)
    for k in range(0,len(db)):
        bestSeq = hit_best_score(db[k],query,m,w)
        if bestSeq != ():
            score = bestSeq[3]
            if score > bestScore or (score == bestScore and bestSeq[2] <res[2]):
                bestScore = score
                res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3],k
    print(bestSeq)
    if bestScore < 0: return ()
    else: return res

def read_database_fasta(filename):
    db = []
    with open(filename, "r") as file:
        for seq_record in SeqIO.parse(file, "fasta"):
            db.append(seq_record)
    return db

def best_alignment_fasta(db,query,w):
    m = build_map(query,w)
    bestScore = -1.0
    res = []#(0,0,0,0,0)
    for k in range(0,len(db)):
        bestSeq = hit_best_score(str(db[k].seq),query,m,w)
        if bestSeq != ():
            score = bestSeq[3]
            if score > bestScore or (score == bestScore and bestSeq[2] <res[2]):
                bestScore = score
                res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3],k
    if bestScore < 0: return ()
    else: return set(res)   

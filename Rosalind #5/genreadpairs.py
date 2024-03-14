#CS321 Generate (k,d)mers from pattern

def pairedcomposition(pattern,k,d):
    kdmers=[]
    for i in range(len(pattern)-(2*k+d)+1):
        readpair = (pattern[i:i+k],pattern[i+k+d:i+k+d+k])
        kdmers.append(readpair)
    kdmers=sorted(kdmers)
    return kdmers


def prefix(pairedcomp):
    k = len(pairedcomp[0][0])
    prefix=dict()
    suffix = dict()
    
    for i in pairedcomp:
        pref = (i[0][:k-1],i[1][:k-1])
        
        suff = (i[0][1:],i[1][1:])
        
        if i in prefix:
            prefix[i].append(pref)
        else:
            prefix[i]=[pref]
        if i in suffix:
            suffix[i].append(suff)
        else:
            suffix[i]=[suff]
    return (prefix)



def findfirstnode(prefix,suffix):
    for i in prefix:
        exists = False
        
        
        for suff in suffix:

            if prefix[i] ==suffix[suff]:
                exists = True
                break
        if (not exists):
            start = i
            break
    return start
        

def paired_dbgraph2(pairedkmers):
    
    k =len(pairedkmers[0][0])
    prefix_list = prefix(pairedkmers)
    dictionary=dict()
    
    
    for i in prefix_list:
        for j in prefix_list[i]:
            lists=[]
            for m in pairedkmers:
                if j[0] == m[0][:k-1] and j[1] == m[1][:k-1]:
                    lists.append([m[0][1:],m[1][1:]])
            
        dictionary[j]=lists
        
    for i in dictionary:
        print(str(i) + "-->" +str(dictionary[i]))
    return dictionary


def findfirst(adjlist):
    for i in adjlist:
        first  = True
        for j in adjlist:
            for k in adjlist[j]:
                if tuple(k)==i:
                    first=False
                    break
        if first ==True:
            start=i
    return start
    
                
                    
            
def reconstruct_from_paireddbgraph(adjlist):
    f = findfirst(adjlist)
    
    
            
            
         
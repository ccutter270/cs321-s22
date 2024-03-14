#CS321 Generate (k,d)mers from pattern

def pairedcomposition(pattern,k,d):
    #given a sequence (pattern) generate the (k,d)mers
    kdmers=[]
    for i in range(len(pattern)-(2*k+d)+1):
        readpair = (pattern[i:i+k],pattern[i+k+d:i+k+d+k])
        kdmers.append(readpair)
    
    
    return kdmers


def prefix(pairedcomp,k):
    #generates the prefixes of each pair
    
  
    prefix=dict()
    
    #the start node is equal to the prefixes of the first pair
    start = (pairedcomp[0][0][:k-1],pairedcomp[0][1][:k-1])
    
    #the end node is the suffixes of the last pair
    end = (pairedcomp[-1][0][1:],pairedcomp[-1][1][1:])
    
    for i in pairedcomp:
        pref = (i[0][:k-1],i[1][:k-1])
        
        if i in prefix:
            prefix[i].append(pref)
        else:
            prefix[i]=[pref]
    
    return prefix,start,end


def paired_dbgraph2(pairedkmers,k):
    #generates the paired debruijn graph from the paired (k,d)mers
    
    k = int(k)
    prefix_list,start,end = prefix(pairedkmers,k)
    dictionary=dict()
    
    
    for i in prefix_list:
        for j in prefix_list[i]:
            lists=[]
            for m in pairedkmers:
                if j[0] == m[0][:k-1] and j[1] == m[1][:k-1]:
                    lists.append([m[0][1:],m[1][1:]])
            
        dictionary[j]=lists
        

    return dictionary,start,end                   
            
def demo():
    #this was our tester function
    k,d,m = readpairs('rosalindsample.txt')
    p = paired_dbgraph2(m,k)
    import julia_eulerian as je
    result = je.create_cycle(p)
    balgraph = result[0]
    start = result[2]
    end = result[1]
    cycle = je.eulerian(balgraph)        
            
    
    



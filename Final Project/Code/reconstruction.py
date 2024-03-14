#Mia Tarantola, Julia Fairbank, Caroline Cutter
#Find an Eulerian Path in a paired debruijn grapah

#importing program that includes functions necessary to generate a paired debruijn graph
import genreadpairs
import random

def readgenome(filename):
    #generates the genome from a txt file

    with open(filename) as file:
        wholegenome=''
        for line in file:
            newline = line.strip('\n')
            wholegenome+=newline
    return(wholegenome)


def create_cycle(graph,start,end):
    """creating the edge from end to start node"""
    
    balance = {}
    in_val = -1
    out_val = 1
    from_val=''
    to_val=''
    
    for val in graph:
        #initializing balance dictionary with values equal to 0
        balance[val] = 0
  
    
    for value in graph:
        #iterating through each node in the graph
        out_edges_list = graph.get(value)
        balance[value] = balance.get(value) + out_val*(len(out_edges_list))
        #filling in graph values 
        for out_edge in out_edges_list:
            balance[tuple(out_edge)] = balance.get(tuple(out_edge), 0) + in_val
    
    
    #checking for unbalanced nodes
    for item in balance:
        if balance.get(item) < 0:
            from_val = item
        elif balance.get(item) > 0:
            to_val = item
    
    #was writing to errors inorder to debug balanced/unbalenced nodes
#     with open('balgraph.txt','w') as file:
#         for i in graph:
#             file.write(str(i)+'\t' + str(balance[i]) + '\n')

    #if the original db graph is cyclic - might not need anymore
    
    if from_val !='':
            
        if from_val in graph:
            graph[from_val].append(to_val)
        else:
            graph[from_val] = [to_val]
    else:
        from_val = end
        to_val = start
        
    
    return(graph, from_val, to_val)

def eulerian(graph):
    'creating eulerian cycle'
    nodes = list(graph.keys())
    cycle = [nodes[0]]
    
    #while there are still nodes left
    while len(graph) > 0:
        
        node = cycle[-1]
        edges = graph.get(tuple(node),[])
        if len(edges) > 0:
            cycle.append(edges.pop())
            if len(edges) == 0:
                graph.pop(tuple(node)) #deleting nodes with no edges
        else:
            for i, new_start in enumerate(cycle):
                if tuple(new_start) in graph:
                    cycle = cycle[i:] + cycle[1:i+1]
                    break
    new_cycle = []
    for i in cycle:
        new_cycle.append(tuple(i))
    return new_cycle
            
def rearrange(cycle, start, end):
    #breaking the cycle and turning into path
    #breaking between start and end node
    for i in range(len(cycle)-1):
        if cycle[i] == end and cycle[(i+1)] == start:
            starting_point = i+1
            ending_point = i
    newcycle = cycle[starting_point::] + cycle[1:ending_point+1]
           
    return newcycle

def getstring(order,k,d):
    #reconstruct sequence from paired db graph
    firsthalf_string = order[0][0]
    secondhalf_string = order[0][1]
    for i in range(1,len(order)):
        #using kmers in position 0
        firsthalf_string+=order[i][0][-1]
        
        #using kmers in position 1
        secondhalf_string +=order[i][1][-1]
  

    if(firsthalf_string[k+d::]!=secondhalf_string[:len(secondhalf_string)-(k+d)]):
     
        print("not a valid path")
    
    
    #checking for overlap
    k = 0
    for i in range(1, len(secondhalf_string)):   
        if firsthalf_string.endswith(secondhalf_string[:i]):
            k = i
   
    #adding part that doesnt overlapa
    merged = firsthalf_string + secondhalf_string[k:]
    return merged

    
    
        
    
    
    
def put_together():
    #TESTER FUNCTION
    #generate paired composition - need to undo hardcode
    m = genreadpairs.pairedcomposition("TAATGCCATGGGATGTT",3,1)
    
    
    #generate paired db graph
    n = genreadpairs.paired_dbgraph2(m)
    for i in n:
        print(str(i)+ "->" +str(n[i]))
    
    #create edge from end to start node
    result = create_cycle(n)
    balanced_graph = result[0]
    
    end = result[1]
    start=result[2]
    
    #create eulerian cycle
    cycle = eulerian(balanced_graph)
    
    #break cycle --> eulerian path
    answer = rearrange(cycle,start,end)
    
    #path --> reconstructed sequence
    print(getstring(answer,k,d))

#sweep k and sweep d compare n50
#randomly select portions of the genome
    
def collect_data(filename):
    whole_seq = readgenome(filename)
    reads=[]
    
    #sample number
    z = 0
    print(str(len(reads))+'num of reads')
#     for i in range(0,166163,1000):
#         reads.append(whole_seq[i:i+1000+1])
        
    #randomly generate 10 reads of length 1000 from the original genome
    for i in range(10):
        startindex = random.randint(0,len(whole_seq)-1001)
        reads.append(whole_seq[startindex:startindex+1000])
        
    #initialize writing to csv file    
    with open("cs321final-data.csv","w") as file:
        
        #excel sheet header
        print('sample','k','d','yes_no (1_0)',file=file,sep=',')
        for read in reads:
            
            for k in range(3,25,2):
                for d in range(1,20,2):
                    #initialize the reconstructed string and read to mismatch
                    yes_no = 0
                    
                    #generate paired composition of sequence
                    pairedcomp = genreadpairs.pairedcomposition(read,k,d)
                    
                    #generate paired db graph
                    dbgraph,start,end = genreadpairs.paired_dbgraph2(pairedcomp,k)
                    
                    #don't anlyze naturally cyclic dbgraphs
                    if (start==end):
                        break
                    
                    #create edge from end to start node
                    result = create_cycle(dbgraph,start,end)
                    balanced_graph = result[0]
                    end = result[1]
                    start=result[2]
                    #create eulerian cycle
                    cycle = eulerian(balanced_graph)
                    #break cycle --> eulerian path
                    answer = rearrange(cycle,start,end)
                    #path --> reconstructed sequence
                    reconstructedstring = getstring(answer,k,d)
                    
                    if reconstructedstring == read:
                        yes_no = 1
                    
                    #write entry to csv file
                    print(str(z),str(k),str(d),str(yes_no),file=file,sep=",")
                    
            z+=1

def demo():
    
    #TESTER FUNCTION
    import genreadpairs as gpr
    k,d,m = gpr.readpairs('rosalindsample.txt')
    k = int(k)
    d = int(d)
    p,start,end = gpr.paired_dbgraph2(m,k)
    
    result = create_cycle(p,start,end)
    balgraph = result[0]
    start = result[2]
    end = result[1]
    
    cycle = eulerian(balgraph)
    answer = rearrange(cycle,start,end)
    return getstring(answer,k,d)
    
    
     
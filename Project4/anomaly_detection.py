from os import listdir
from statistics import median
import matplotlib.pyplot as plt
import networkx as nx
import hashlib
import sys

def main():
    file_list = get_file_list()
    simHash_signature_list = get_simhash_signature_list( file_list )    
    sim_list = []
    for x,y in zip( simHash_signature_list, simHash_signature_list[1:] ):
        sim = calculate_similarity(x,y)
        sim_list.append(sim)

    mr_multiplier = 3
    if len(sys.argv) == 4:
        mr_multiplier = float(sys.argv[3])

    threshold = calculateThreshold(sim_list, mr_multiplier)
    out_put_file = sys.argv[2]
    plot_result(sim_list, threshold, out_put_file+".png" ,mr_multiplier)
    anomalous_point = find_anomalous_points(sim_list, threshold)
    write_output(anomalous_point, out_put_file)
    
def get_file_list():
    if len(sys.argv) < 3:
        print("Follow the running instruction:")
        exit()

    data_directory = sys.argv[1] + "/"
    #Sort the file list based on integer before underscore
    file_list = sorted( listdir(data_directory), key=lambda item: (int(item.partition('_')[0])) )
    return [ data_directory+x for x in file_list ]

def get_simhash_signature_list(file_list):
    simhash_list = []
    G = nx.DiGraph()
    
    for file_path in file_list:
        graph_file = open(file_path,'r')
        G.clear()
        for line in graph_file:
            n1,n2 = [int(x) for x in line.split(' ')]
            G.add_edge(n1,n2)
        graph_file.close()
        
        feature_set = create_graph_feature_set(G)

        sim_hash_vector = sim_hash( feature_set )
        simhash_list.append( sim_hash_vector )
    return simhash_list

def create_graph_feature_set(G):  
    page_rank =  nx.pagerank(G) 
    feature_set = [ (str(ti),wi) for ti,wi in page_rank.items() ]

    for edge in G.edges():
        ti = str(edge[0]) + ' ' + str(edge[1])
        out_degree = G.out_degree(edge[0])       
        wi = page_rank[ edge[0] ] /out_degree 
        feature_set.append( (ti,wi) )
    return feature_set

def sim_hash( feature_set ):
    V = [0]*128
    for ti, wi in feature_set:
        hex_dig = hashlib.md5( ti.encode('utf-8') ).hexdigest()
        bin_repr = bin( int( hex_dig, 16) )[2:].zfill(128)       
        for i in range(128):
            if bin_repr[i] == '1':
                V[i] = V[i] + wi  
            else: 
                V[i] = V[i] - wi

    for i in range(128):
        if V[i] > 0: V[i] = 1  
        else: V[i] = 0
    return(V)

def calculate_similarity(x,y):
    assert len(x) == len(y)
    b = len(x)
    hamming_dist = 0
    for i in range(b):
        if x[i] != y[i]: hamming_dist += 1
    return(1 - hamming_dist/b) 

def calculateThreshold(similarity_list, multiple_of_mr):
    m = median(similarity_list)
    n = len(similarity_list) 
    if n < 2: 
        return(m)
    mr_sum = 0
    for i in range(1,n):
        mr_i = abs(similarity_list[i] - similarity_list[i-1]) # mr(i) = | x(i+1) - x(i) |
        mr_sum += mr_i  
    mr = mr_sum/(n-1)
    return( {"lower": m - multiple_of_mr * mr, "median": m, "mr" : mr} )

def plot_result(sim_list, threshold, filename, mr_multiplier):
    """A simple function to plot the the points and the threshold lines"""
    plt.plot(range( len(sim_list) ), sim_list, 'ro', color = '0.75' )
    plt.title('Threshold: (Median - ' + str(mr_multiplier) + ' * MR)' )
    plt.xlabel('Index')
    plt.ylabel('Similarity')
    plt.grid(True)   
    line1, = plt.plot([0, len(sim_list)], [threshold["lower"], threshold["lower"]], 'k--', lw=1, alpha=0.75)
    plt.legend([line1], ['Threshold'])
    plt.savefig(filename)

def find_anomalous_points(sim_list, threshold):
    anomalous_point = []
    for index, similarity in enumerate( sim_list[:-1] ): 
        next_similarity = sim_list[index+1]
        if ( (similarity < threshold["lower"] and next_similarity < threshold["lower"]) ):
            distance_from_median = abs( similarity - threshold["median"] )
            anomalous_point.append( (index+1, distance_from_median) )
    return(anomalous_point)

def write_output(anomalous_points, out_put_file):
    number_of_anomalies = len(anomalous_points)
    out_file = open(out_put_file + ".txt" , 'w+') 
    out_file.write(str(number_of_anomalies))

    anomalous_points.sort(key= lambda item: item[1], reverse=True)

    for i in range(number_of_anomalies):
        out_file.write( '\n' + str( anomalous_points[i][0] ) )

    out_file.close()

main()
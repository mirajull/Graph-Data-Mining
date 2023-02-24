from igraph import *
from os.path import join
from os import listdir
import scipy.linalg as sc
from pylab import *
import numpy as np
import random
import time

def highest_eigen_value(g):
    length = len(g.vs)
    eg_val = sc.eigh(g.get_adjacency().data,eigvals_only=True,eigvals=(length-1,length-1))
    return eg_val[0]

def shield_score(g):
    length = len(g.vs)
    A = np.matrix(g.get_adjacency().data)
    eg = sc.eigh(A,eigvals_only=False,eigvals=(length-1,length-1))
    lamda = eg[0][0]
    u = eg[1]
    S=[]
    v=[0,0,0,0,0,0,0,0,0,0]
    score = [0,0,0,0,0,0,0,0,0,0]
    for j in range(length):
        v[j] = (2*lamda - A[j,j])*u[j][0]*u[j][0]
    for i in range(3):
        B = A[:,S]
        b = np.dot(B,u[S])
        for j in range(length):
            if j in S:
                score[j] = -1
            else:
                score[j] = v[j] - 2 *b[j]*u[j]
        l = np.argmax(score)
        S.append(l)
    return S


def effective_strength(high_eig,beta,delta):
    return round((beta/delta)*high_eig , 3)

def plot_beta(high_eig,delta,filename,save):
    betas = np.arange(0.,1.,0.1)
    strengths = [effective_strength(high_eig,b,delta) for b in betas]
    min_beta = delta/high_eig

    figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(betas,strengths, "-o")
    axhline(y=min_beta, ls='-', c='r',
            label='Minimum Transimission Probability: %.2f'%(min_beta),
            lw=2)
    plt.grid(True)
    plt.legend(loc='best')
    plt.title("Virus Propogation")
    plt.xlabel('Beta')
    plt.ylabel('Effective Strength')
    if save:
        savefig(join('graph', filename+"_beta"+str(time.time())+".png"),bbox_inches='tight')
    else:
        plt.show()
    print("Minimum transmission probability: {}".format(min_beta))


def plot_delta(high_eig,beta,filename,save):
    deltas = np.arange(0.,1.,0.1)
    strengths = [effective_strength(high_eig,beta,d) for d in deltas]
    max_delta = beta*high_eig

    #plot the results
    figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(deltas,strengths,'o-')
    axhline(y=max_delta, ls='-', c='r',
            label='Maximum Healing Probability: %.3f'%(max_delta),
            lw=2)
    plt.grid(True)
    plt.legend(loc='best')
    plt.title("Virus Propogation")
    plt.xlabel('Delta')
    plt.ylabel('Effective Strength')
    if save:
        savefig(join('graph', filename+"_delta"+str(time.time())+".png"),bbox_inches='tight')
    else:
        plt.show()
    print("Maximum healing probability: {}".format(max_delta))


def get_infected_neighbors(g,n,beta):
    infected =  [k for k in g.neighbors(g.vs[n]) if np.random.uniform(0,1) <= beta]
    return infected


def simulate(g ,beta,delta,t,c):
    samples = random.sample(range(g.vs[0].index,g.vs[len(g.vs)-1].index),c)
    infected = samples
    num_infected = [len(samples)]
    for i in range(1, t):
        infected_at_t = sum([get_infected_neighbors(g,n,beta) for n in infected],[])
        cured =  [k for k in infected if np.random.uniform(0,1) <= delta]
        infected = list(set(infected)-set(cured))
        infected = list(set(infected + infected_at_t))
        num_infected.append(len(infected))
    return num_infected

def run_simulation(g,beta,delta,t,c,runs):
    simulations=[]
    for i in range(0,runs):
        simulations = simulations + [simulate(g,beta,delta,t,c)]
    frac_affected = []
    for i in range(0,t):
        frac_affected.append(mean([simulations[j][i] for j in range(runs)])/len(g.vs))
    return frac_affected

def random_immunize(g,k):
    g1 = g.copy()
    nodes = random.sample(range(len(g1.vs)),k)
    g1.delete_vertices(nodes)
    return g1

def high_degree_immunize(g,k):
    g1 = g.copy()
    degree = {v: g1.degree(v) for v in g1.vs}
    sorted_deg = sorted(degree.keys(),key = lambda s : degree[s],reverse = True)
    g1.delete_vertices(sorted_deg[:k])
    return g1


def high_degree_immunize_iteratively(g,k):
    g1 = g.copy()
    for i in range(k):
        g1 = high_degree_immunize(g1,1)
    return g1


def largest_eig_vec_immunize(g,k):
    g1 = g.copy()
    largest_eig = sc.eigh(g1.get_adjacency().data,eigvals=(len(g1.vs)-1,len(g1.vs)-1))
    eig_vec = {i: largest_eig[1][i] for i in range(len(largest_eig[1]))}
    sorted_deg = sorted(eig_vec.keys(),key = lambda s : abs(eig_vec[s][0]),reverse = True)
    max_eig_vertices =  [g1.vs[i] for i in sorted_deg][:k]
    g1.delete_vertices(max_eig_vertices)
    return g1



def plot_simulation(f_infected,t,filename,save):
    t = range(0,t)
    figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(t,f_infected, "-o")
    plt.grid(True)
    plt.title("Virus Propogation")
    plt.xlabel('Time')
    plt.ylabel('Fraction of Infected nodes')
    if save:
        savefig(join('graph', filename+"_virus_simulation"+str(time.time())+".png"),bbox_inches='tight')
    else:
        plt.show()


def plot_k_variations(g,k_list,beta,delta,policy,filename,save):
    g1 = g.copy()
    eff_strength = []
    for k in k_list:
        if policy == "A":
            g1 = random_immunize(g1,k)
        if policy == "B":
            g1 = high_degree_immunize(g1,k)
        if policy == "C":
            g1 = high_degree_immunize_iteratively(g1,k)
        if policy == "D":
            g1 = largest_eig_vec_immunize(g1,k)
        h_eigen = highest_eigen_value(g1)
        eff_strength.append(effective_strength(h_eigen,beta,delta))
        g1 = g.copy()
    figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(k_list,eff_strength, "-o")
    plt.grid(True)
    plt.title("Immunization chart and Effective Strength")
    plt.xlabel('Number of vaccines')
    plt.ylabel('Effective Strength')
    if save:
        savefig(join('graph', filename+"_immunization_simulation"+str(time.time())+".png"),bbox_inches='tight')
    else:
        plt.show()


def get_system_matrix(graphs,beta,delta):
    S=[]
    for graph in graphs:
        Si = (1-delta)*np.identity(len(graph.vs)) + beta*np.array(graph.get_adjacency().data)
        S.append(Si)
    system_matrix = S[0]
    for i in range(1,len(S)):
        system_matrix = system_matrix.dot(S[i])
    return system_matrix


def get_graph(file):
    g = Graph()
    with open(file) as fi:
            v,e = fi.next().split()
            e_list = [(int(line.split()[0]) , int(line.split()[1])) for line in list(fi)]
            g.add_vertices(int(v))
            g.add_edges(e_list)
    return g


if __name__ == '__main__':
    beta1 = 0.20
    beta2 = 0.01
    delta1 = 0.70
    delta2 = 0.60
    k1 = 200

    g = get_graph(sys.argv[1])
    highest_e_val = highest_eigen_value(g)
    eff_strength1 = effective_strength(highest_e_val,beta1,delta1)
    print("The effective strength of the virus for Beta = {} and Delta = {} is {}".format(beta1,delta1,eff_strength1))
    if eff_strength1 > 1:
        print("The infection will spread across the network.")
    else:
        print("The virus will die quickly.")

    eff_strength2 = effective_strength(highest_e_val,beta2,delta2)
    print("The effective strength of the virus for Beta = {} and Delta = {} is {} ".format(beta2,delta2,eff_strength2))
    if eff_strength2 > 1:
        print("The infection will spread across the network.")
    else:
        print("The virus will die quickly.")

    # plot_beta(highest_e_val,delta1,out_file,True)
    # plot_delta(highest_e_val,beta1,out_file,True)
    # plot_beta(highest_e_val,delta1,out_file,True)
    # plot_delta(highest_e_val,beta2,out_file,True)

    # frac_infected = run_simulation(g,beta1,delta1,100,len(g.vs)/10,10)
    # plot_simulation(frac_infected,100,out_file,True)
    # frac_infected = run_simulation(g,beta2,delta2,100,len(g.vs)/10,10)
    # plot_simulation(frac_infected,100,out_file,True)

#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import pickle
import sys
import textwrap
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, "r") as file_in:
        line = file_in.readline()
        while line != "":
            line = file_in.readline()
            yield line.strip()
            file_in.readline()
            file_in.readline()
            file_in.readline()

def cut_kmer(read, kmer_size):
    nb_iter = len(read) - (kmer_size-1)
    for i in range(nb_iter):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    read = "".join(list(read_fastq(fastq_file)))
    list_kmer = cut_kmer(read, kmer_size)
    kmer_dict = {}

    for kmer in list_kmer:
        kmer_dict[kmer] = kmer_dict.get(kmer,0)+1

    return kmer_dict

def build_graph(kmer_dict):
    digraph = nx.DiGraph()
    for kmer, weight in kmer_dict.items():
        digraph.add_edge(kmer[:-1], kmer[1:], weight=weight)
    return digraph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)
        elif delete_entry_node and delete_sink_node == False:
            graph.remove_nodes_from(path[:-1])
        elif delete_entry_node == False and delete_sink_node:
            graph.remove_nodes_from(path[1:])
        else:
            graph.remove_nodes_from(path[1:-1])
    return graph

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    std_avg_weight = statistics.stdev(weight_avg_list)
    if std_avg_weight > 0 :
        best_path = weight_avg_list.index(max(weight_avg_list))
    else :
        std_len = statistics.stdev(path_length)
        if std_len > 0:
            best_path = path_length.index(max(path_length))
        else:
            best_path = randint(len(path_list)-1)
    paths_to_remove = path_list[:best_path]+path_list[best_path+1:]
    remove_paths(graph,paths_to_remove,delete_entry_node,delete_sink_node)
    return graph        

def path_average_weight(graph, path):
    """Compute the weight of a path"""
    return statistics.mean([d["weight"] for (u,v,d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = list(nx.all_simple_paths(graph, source = ancestor_node, target = descendant_node))

    weight_avg_list = []
    path_length = []

    for path in path_list:
        weight_avg_list.append(path_average_weight(graph, path))
        path_length.append(len(path))

    return select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False)

def simplify_bubbles(graph):
    bubble = False 
    for node in graph:
        list_preds = list(graph.predecessors(node))
        if len(list_preds) > 1:
            for i,first_pred in enumerate(list_preds):
                list_second_preds = list_preds[:i]+list_preds[i+1:]
                for second_pred in list_second_preds:
                    ancestor_node = nx.lowest_common_ancestor(graph, first_pred, second_pred)
                    if ancestor_node:
                        bubble = True
                        break
        if bubble == True:
            break
   # La simplification ayant pour conséquence de supprimer des noeuds du hash
   # Une approche récursive est nécessaire avec networkx
    if bubble:
        graph = simplify_bubbles(solve_bubble(graph,ancestor_node, node))

    return graph

def solve_entry_tips(graph, starting_nodes):
    path_list = []
    pred_start = 0 
    for node in graph:
        list_preds = list(graph.predecessors(node))
        if len(list_preds) > 1:
            for start_node in starting_nodes:
                path = list(nx.all_simple_paths(graph,source=start_node,target=node))[0]
                if path:
                    pred_start += 1
                    path_list.append(path)
                    if pred_start == 2 :
                        break
        if pred_start == 2:
            break
    
    if pred_start == 2 :
        weight_avg_list = []
        path_length = []

        for path in path_list:
            print(path)
            weight_avg_list.append(path_average_weight(graph, path))
            path_length.append(len(path))

        graph = select_best_path(graph, path_list, path_length, weight_avg_list, 
                        delete_entry_node=True, delete_sink_node=False)
        ending_nodes = get_starting_nodes(graph)
        graph = solve_entry_tips(graph,ending_nodes)
    return graph

def solve_out_tips(graph, ending_nodes):
    path_list = []
    for i,first_node in enumerate(ending_nodes):
        other_ending_nodes = ending_nodes[:i]+ending_nodes[i+1:]
        for second_node in other_ending_nodes:
            ancestor_node = nx.lowest_common_ancestor(graph,first_node,second_node)
            if ancestor_node:
                first_path = list(nx.all_simple_paths(graph,source=ancestor_node,target=first_node))[0]
                path_list.append(first_path)
                second_path = list(nx.all_simple_paths(graph,source=ancestor_node,target=second_node))[0]
                path_list.append(second_path)
                break
        if path_list:
            break
    
    if path_list:
        weight_avg_list = []
        path_length = []

        for path in path_list:
            weight_avg_list.append(path_average_weight(graph, path))
            path_length.append(len(path))

        graph = select_best_path(graph, path_list, path_length, weight_avg_list, 
                        delete_entry_node=False, delete_sink_node=True)
        ending_nodes = get_sink_nodes(graph)
        graph = solve_out_tips(graph,ending_nodes)
    return graph

def get_starting_nodes(graph):
    starting_nodes =[]
    for node in graph:
        pred = list(graph.predecessors(node))
        if not(pred):
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(graph):
    sink_nodes =[]
    for node in graph:
        pred = list(graph.successors(node))
        if not(pred):
            sink_nodes.append(node)
    return sink_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for start in starting_nodes:
        for end in ending_nodes:
            if nx.has_path(graph, start, end):
                seq = ""
                for path in nx.all_simple_paths(graph, start, end):
                    seq += path[0]
                    for node in path[1:]:
                        seq += node[-1]
                    contigs.append((seq, len(seq)))
    return contigs

def save_contigs(contigs_list, output_file):
    with open(output_file,"wt") as file:
        for i,contig in enumerate(contigs_list):
            seq = contig[0]
            len = contig[1]
            file.write(f">contig_{i} len={len}\n{textwrap.fill(seq)}\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()

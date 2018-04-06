from __future__ import print_function, division

from csv import reader as csv_reader
from itertools import chain
import networkx as nx
import phylonetwork as ph

import argparse


def read_enewick(s):
    return ph.PhyloNetwork(eNewick=s)

def iter_length(g):
    c = 0
    for _ in g:
        c += 1
    return c

def nodes_below(net, u):
    yield u

    u_successors = nx.dfs_successors(net, u)
    nodes_below_u = chain.from_iterable(u_successors.values())
    for v in nodes_below_u:
        yield v

def kappa(net, u):
    if net.is_leaf(u):
        return 1

    return iter_length(filter(net.is_leaf, nodes_below(net, u)))

def cache_kappa(net):
    values = {u:kappa(net, u) for u in net}
    return values.__getitem__

def fair_proportion(net, u, weight='weight', kappa_fn=None):
    if kappa_fn is None:
        kappa_fn = cache_kappa(net)

    if isinstance(weight, str):
        weight_key = weight
        weight = lambda n, w, v: n.edge[w][v][weight_key]

    rev_net = nx.reverse(net)
    nodes_above_u = nodes_below(rev_net, u)
    edges_above_u = ((w, v) for v in nodes_above_u for w in rev_net[v])
    return sum(weight(net, w, v) / kappa_fn(v) for w, v in edges_above_u)



def read_weights_file(f):
    return {(u, v): float(w) for u,v,w in csv_reader(f, skipinitialspace=True)}

FP_USAGE_EXAMPLE = u"""
Usage example:

  # Computing the fair proportion of nodes 1, 2, 3, 4 and 5
  $ %(prog)s fair-proportion '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 1 2 3 4 5 <<EOF
      H,4,0.5
      a,c,0.3
      b,d,0.4
      r,a,0.1
      r,b,0.2
  EOF
  1\t0.0333333333333
  2\t0.0666666666667
  3\t0.183333333333
  4\t0.95
  5\t0.266666666667
\u00a0
"""

def get_args_network(args):
    net = args.network
    weights = read_weights_file(args.weights_file)

    for u in net:
        for v in net[u]:
            labs = net.label(u), net.label(v)
            if labs in weights:
                net.edge[u][v]['weight'] = weights[labs]
            else:
                net.edge[u][v]['weight'] = 0

    return net

def main_fair_proportion(args):
    net = get_args_network(args)
    kappa_fn = cache_kappa(net)

    for node_name in args.nodes:
        print('{}\t{}'.format(node_name, fair_proportion(net, net.node_by_taxa(node_name), kappa_fn=kappa_fn)))

def main():
    def phylonetwork_argument(s):
        try:
            return read_enewick(s)
        except ph.classes.MalformedNewickException as e:
            raise argparse.ArgumentTypeError(str(e))

    def add_network_arg(parser):
        parser.add_argument(
            'network',
            help='Phylogenetic network in eNewick format (weights not supported)',
            type=phylonetwork_argument)

    def add_weights_file_arg(parser):
        parser.add_argument(
            'weights_file',
            help='Edge weights file. It should be in CSV format with three columns and no header: node1, node2, weight.',
            type=argparse.FileType('r'))

    parser = argparse.ArgumentParser(
        description='Compute phylogenetic network indices')

    subparsers = parser.add_subparsers(title='actions')

    fp_parser = subparsers.add_parser('fair-proportion',
        help='Compute the fair proportion of a node',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=FP_USAGE_EXAMPLE)

    fp_parser.set_defaults(main_fn=main_fair_proportion)

    add_network_arg(fp_parser)
    add_weights_file_arg(fp_parser)
    fp_parser.add_argument('nodes', help='Target node', type=str, nargs='+')

    args = parser.parse_args()
    args.main_fn(args)


if __name__ == '__main__':
    main()

    #net = read_enewick('((1,(3,#H1)c:0.3)a:0.1,(2,((4:0.5)H#H1,5)d:0.4)b:0.2)r;')
    #weights = {('H','4'): 0.5, ('a','c'): 0.3, ('b','d'): 0.4, ('r','a'): 0.1, ('r','b'): 0.2}
    # H,4,0.5
    # a,c,0.3
    # b,d,0.4
    # r,a,0.1
    # r,b,0.2

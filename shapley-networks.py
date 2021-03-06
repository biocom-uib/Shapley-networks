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

def default_weight_fn(net, w, v):
    return net.edge[w][v].get('length', 0)

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

def cophenetic_value(net, us, weight_fn=default_weight_fn):
    if not us:
        return 0

    rev_net = nx.reverse(net)
    nodes_above_all_us = frozenset.intersection(*[frozenset(nodes_below(rev_net, u)) for u in us])
    edges_above_all_us = ((w, v) for v in nodes_above_all_us for w in rev_net[v])
    return sum(weight_fn(net, w, v) for w, v in edges_above_all_us)

def rooted_phylogenetic_diversity(net, us, weight_fn=default_weight_fn):
    rev_net = nx.reverse(net)
    nodes_above_us = frozenset().union(*[frozenset(nodes_below(rev_net, u)) for u in us])
    edges_above_us = ((w, v) for v in nodes_above_us for w in rev_net[v])
    return sum(weight_fn(net, w, v) for w, v in edges_above_us)

def unrooted_phylogenetic_diversity(net, us, weight_fn=default_weight_fn):
    return rooted_phylogenetic_diversity(net, us, weight_fn) - cophenetic_value(net, us, weight_fn)

def fair_proportion(net, u, weight_fn=default_weight_fn, kappa_fn=None):
    if kappa_fn is None:
        kappa_fn = cache_kappa(net)

    rev_net = nx.reverse(net)
    nodes_above_u = nodes_below(rev_net, u)
    edges_above_u = ((w, v) for v in nodes_above_u for w in rev_net[v])
    return sum(weight_fn(net, w, v) / kappa_fn(v) for w, v in edges_above_u)

def cophenetic_shapley_value(net, u, weight_fn=default_weight_fn, kappa_fn=None):
    if kappa_fn is None:
        kappa_fn = cache_kappa(net)

    leaves = net.leaves()
    n = len(leaves)
    rev_net = nx.reverse(net)
    nodes_above_u = nodes_below(rev_net, u)
    edges_above_u = frozenset((w, v) for v in nodes_above_u for w in rev_net[v])
    edges_not_above_u = (e for e in net.edges() if e not in edges_above_u)

    total_weight = sum(weight_fn(net, *e) for e in net.edges())

    return total_weight/n - sum(weight_fn(net, w, v) / (n - kappa_fn(v))
                                  for w, v in edges_not_above_u)

def unrooted_shapley_value(net, u, weight_fn=default_weight_fn, kappa_fn=None):
    if kappa_fn is None:
        kappa_fn = cache_kappa(net)

    sv = fair_proportion(net, u, weight_fn=weight_fn, kappa_fn=kappa_fn)
    coph_sv = cophenetic_shapley_value(net, u, weight_fn=weight_fn, kappa_fn=kappa_fn)
    return sv - coph_sv

COPH_USAGE_EXAMPLE = u"""
Usage example:

  # Computing the cophenetic value of {3,4}
  $ %(prog)s cophenetic-value '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 1 2 3 <<EOF
      c,H,0.5
      H,4,0.5
      a,c,0.3
      b,d,0.4
      r,a,0.1
      r,b,0.2
  EOF
  0.4
\u00a0
"""

rPSD_USAGE_EXAMPLE = u"""
Usage example:

  # Computing the rooted phylogenetic diversity of {3,4}
  $ %(prog)s rpsd '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 3 4 <<EOF
      H,4,0.5
      a,c,0.3
      b,d,0.4
      r,a,0.1
      r,b,0.2
  EOF
  1.5
\u00a0
"""

uPSD_USAGE_EXAMPLE = u"""
Usage example:

  # Computing the unrooted phylogenetic diversity of nodes {3,4}
  $ %(prog)s upsd '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 3 4 <<EOF
      H,4,0.5
      a,c,0.3
      b,d,0.4
      r,a,0.1
      r,b,0.2
  EOF
  1.1
\u00a0
"""

FP_USAGE_EXAMPLE = u"""
Usage example:

  # Computing the fair proportion of all the leaves
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

COPH_SV_USAGE_EXAMPLE = u"""
Usage example:

  # Computing the cophenetic shapley value of all the leaves
  $ %(prog)s cophenetic-shapley-value '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 1 2 3 4 5 <<EOF
      H,4,0.5
      a,c,0.3
      b,d,0.4
      r,a,0.1
      r,b,0.2
  EOF
  1\t-0.158333333333
  2\t-0.108333333333
  3\t-0.0583333333333
  4\t0.3
  5\t0.025
\u00a0
"""

U_SV_USAGE_EXAMPLE = u"""
Usage example:

  # Computing the unrooted shapley value of all the leaves
  $ %(prog)s unrooted-shapley-value '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 1 2 3 4 5 <<EOF
      H,4,0.5
      a,c,0.3
      b,d,0.4
      r,a,0.1
      r,b,0.2
  EOF
  1\t0.191666666667
  2\t0.175
  3\t0.241666666667
  4\t0.65
  5\t0.241666666667
\u00a0
"""

def main_cophenetic_value(args):
    net = args.network

    print(cophenetic_value(net, [net.node_by_taxa(node_name) for node_name in args.nodes]))

def main_rooted_phylogenetic_diversity(args):
    net = args.network

    print(rooted_phylogenetic_diversity(net, [net.node_by_taxa(node_name) for node_name in args.nodes]))

def main_unrooted_phylogenetic_diversity(args):
    net = args.network

    print(unrooted_phylogenetic_diversity(net, [net.node_by_taxa(node_name) for node_name in args.nodes]))

def main_fair_proportion(args):
    net = read_enewick(args.network.eNewick())
    kappa_fn = cache_kappa(net)

    for node_name in args.nodes:
        print('{}\t{}'.format(node_name, fair_proportion(net, net.node_by_taxa(node_name), kappa_fn=kappa_fn)))

def main_cophenetic_shapley_value(args):
    net = args.network
    kappa_fn = cache_kappa(net)

    for node_name in args.nodes:
        print('{}\t{}'.format(node_name, cophenetic_shapley_value(net, net.node_by_taxa(node_name), kappa_fn=kappa_fn)))

def main_unrooted_shapley_value(args):
    net = args.network
    kappa_fn = cache_kappa(net)

    for node_name in args.nodes:
        print('{}\t{}'.format(node_name, unrooted_shapley_value(net, net.node_by_taxa(node_name), kappa_fn=kappa_fn)))


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

    parser = argparse.ArgumentParser(
        description='Compute phylogenetic network indices')

    subparsers = parser.add_subparsers(title='actions')


    coph_parser = subparsers.add_parser('cophenetic-value',
        help='Compute the cophenetic value of a subnet',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=COPH_USAGE_EXAMPLE)

    add_network_arg(coph_parser)
    coph_parser.add_argument('nodes', help='Bottom nodes in the subnet', type=str, nargs='*')
    coph_parser.set_defaults(main_fn=main_cophenetic_value)

    rpsd_parser = subparsers.add_parser('rpsd',
        help='Compute the (rooted) phylogenetic subnet diversity',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=rPSD_USAGE_EXAMPLE)

    add_network_arg(rpsd_parser)
    rpsd_parser.add_argument('nodes', help='Bottom nodes in the subnet', type=str, nargs='*')
    rpsd_parser.set_defaults(main_fn=main_rooted_phylogenetic_diversity)

    upsd_parser = subparsers.add_parser('upsd',
        help='Compute the (unrooted) phylogenetic subnet diversity',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=uPSD_USAGE_EXAMPLE)

    add_network_arg(upsd_parser)
    upsd_parser.add_argument('nodes', help='Bottom nodes in the subnet', type=str, nargs='*')
    upsd_parser.set_defaults(main_fn=main_unrooted_phylogenetic_diversity)

    fp_parser = subparsers.add_parser('fair-proportion',
        help='Compute the fair proportion of a node',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=FP_USAGE_EXAMPLE)

    add_network_arg(fp_parser)
    fp_parser.add_argument('nodes', help='Target nodes', type=str, nargs='+')
    fp_parser.set_defaults(main_fn=main_fair_proportion)

    coph_sv_parser = subparsers.add_parser('cophenetic-shapley-value',
        help='Compute the cophenetic shapley value of a node',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=COPH_SV_USAGE_EXAMPLE)

    add_network_arg(coph_sv_parser)
    coph_sv_parser.add_argument('nodes', help='Target nodes', type=str, nargs='+')
    coph_sv_parser.set_defaults(main_fn=main_cophenetic_shapley_value)

    u_sv_parser = subparsers.add_parser('unrooted-shapley-value',
        help='Compute the unrooted shapley value of a node',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=U_SV_USAGE_EXAMPLE)

    add_network_arg(u_sv_parser)
    u_sv_parser.add_argument('nodes', help='Target nodes', type=str, nargs='+')
    u_sv_parser.set_defaults(main_fn=main_unrooted_shapley_value)

    args = parser.parse_args()
    args.main_fn(args)


if __name__ == '__main__':
    main()
    #net = read_enewick('((1,(3,#H1)c:0.3)a:0.1,(2,((4:0.5)H#H1,5)d:0.4)b:0.2)r;')

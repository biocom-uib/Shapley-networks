# Shapley networks

This script can be used to compute some indices of a phylogenetic network, such
as the fair proportion of a node or its unrooted shapley values.

Note that all networks are introduced in the [Extended Newick](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-532) format.

## Dependencies

Before using this script, you need to install [PhyloNetwork](https://github.com/bielcardona/PhyloNetworks) `>= 1.2`:
```
pip3 install --user phylonetwork
```

## Computing the phylogenetic subnet diversity

Given the network `((1,(3,#H1)c:0.3)a:0.1,(2,((4:0.5)H#H1,5)d:0.4)b:0.2)r;`,
one can obtain both the rooted and the unrooted phylogenetic subnet diversity
by executing (for instance, with `X = {3,4}`)
```
$ python3 shapley-networks.py rpsd '((1,(3,#H1)c:0.3)a:0.1,(2,((4:0.5)H#H1,5)d:0.4)b:0.2)r;' 3 4
1.5

$ python3 shapley-networks.py upsd '((1,(3,#H1)c:0.3)a:0.1,(2,((4:0.5)H#H1,5)d:0.4)b:0.2)r;' 3 4
1.1
```

Similarly, the cophenetic value is obtained via the `cophenetic-value` command
```
$ python3 shapley-networks.py cophenetic-value '((1,(3,#H1)c:0.3)a:0.1,(2,((4:0.5)H#H1,5)d:0.4)b:0.2)r;' 3 4
0.4
```

## Computing the Shapley values of a node

Given the network `((1,(3,#H1)c:0.3)a:0.1,(2,((4:0.5)H#H1,5)d:0.4)b:0.2)r;`,
the FP of nodes 1, 2, 3, 4 and 5 can be computed as follows:
```
$ python3 shapley-networks.py fair-proportion '((1,(3,#H1)c:0.3)a:0.1,(2,((4:0.5)H#H1,5)d:0.4)b:0.2)r;' 1 2 3 4 5
1       0.03333333333333333
2       0.06666666666666667
3       0.18333333333333332
4       0.95
5       0.26666666666666666
```

Cophenetic Shapley values and the unrooted Shapley values are calculated analogously:
```
$ python3 shapley-networks.py cophenetic-shapley-value '((1,(3,#H1)c:0.3)a:0.1,(2,((4:0.5)H#H1,5)d:0.4)b:0.2)r;' 1 2 3 4 5
1       -0.15833333333333338
2       -0.10833333333333334
3       -0.05833333333333335
4       0.3
5       0.024999999999999967

$ python3 shapley-networks.py unrooted-shapley-value '((1,(3,#H1)c:0.3)a:0.1,(2,((4:0.5)H#H1,5)d:0.4)b:0.2)r;' 1 2 3 4 5
1       0.1916666666666667
2       0.175
3       0.24166666666666667
4       0.6500000000000001
5       0.2416666666666667
```

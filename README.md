# Shapley networks

## Dependencies

Before using this script, you need to install [PhyloNetwork](https://github.com/bielcardona/PhyloNetworks):
```
pip2 install --user phylonetwork
```

## Computing the phylogenetic subnet diversity

Given the network `((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r` with the weights given
below, one can obtain both the rooted and the unrooted phylogenetic subnet
diversity by executing (for instance, with `X = {3,4}`)
```
$ python2 shapley-networks.py rpsd '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 3 4 <<EOF
H,4,0.5
a,c,0.3
b,d,0.4
r,a,0.1
r,b,0.2
EOF
1.5

$ python2 shapley-networks.py upsd '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 3 4 <<EOF
H,4,0.5
a,c,0.3
b,d,0.4
r,a,0.1
r,b,0.2
EOF
1.1
```

Similarly, the cophenetic value is obtained via the `cophenetic-value` command
```
$ python2 shapley-networks.py cophenetic-value '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 3 4 <<EOF
H,4,0.5
a,c,0.3
b,d,0.4
r,a,0.1
r,b,0.2
EOF
0.4
```

## Computing the Shapley values of a node

Given the network `((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r`, the FP of the nodes 1,
2, 3, 4 and 5 can be computed as follows:
```
$ python2 shapley-networks.py fair-proportion '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 1 2 3 4 5 <<EOF
H,4,0.5
a,c,0.3
b,d,0.4
r,a,0.1
r,b,0.2
EOF
1       0.0333333333333
2       0.0666666666667
3       0.183333333333
4       0.95
5       0.266666666667
```

The cophenetic Shapley value and the unrooted Shapley value are calculated analogously:
```
$ python2 shapley-networks.py cophenetic-shapley-value '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 1 2 3 4 5 <<EOF
H,4,0.5
a,c,0.3
b,d,0.4
r,a,0.1
r,b,0.2
EOF
1       -0.158333333333
2       -0.108333333333
3       -0.0583333333333
4       0.3
5       0.025

$ python2 shapley-networks.py unrooted-shapley-value '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 1 2 3 4 5 <<EOF
H,4,0.5
a,c,0.3
b,d,0.4
r,a,0.1
r,b,0.2
EOF
1       0.191666666667
2       0.175
3       0.241666666667
4       0.65
5       0.241666666667
```

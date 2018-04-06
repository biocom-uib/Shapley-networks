# Shapley networks

## Dependencies

`phylonetwork`

## Computing the FP of a network

Given the network ((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r, the FP of the nodes 1, 2,
3, 4 and 5 can be computed as follows:
```
$ python2 shapley-networks.py fair-proportion '((1,(3,#H1)c)a,(2,((4)H#H1,5)d)b)r;' - 1 2 3 4 5 <<EOF
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
```

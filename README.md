# NewGenomeSeqDetection

These codes aim to address the new genome sequence detection problem.

## Getting Started

Some packages and files are needed to perform RAP and RAPCOS. Please fully read this document before implementation, in particular, please pay attention to the **Notice** part.

### Prerequisites

The python packages sys, random, numpy and pandas are required to perform RAP and RAPCOS.

The files ```SeqVer.fasta```, ```integer.xlsx```, ```vertex.xlsx``` are needed for demo codes.

### Examples for RAP

```
python RAP.py
```

where the ```main()``` function can be performed directly. 

```RAP.py``` is used to perform RAP with loss function measuring the distance between natural vector and target natural vector. 

### Examples for RAPCOS

```
python RAPCOS.py
```

where the ```main()``` function can be performed directly. 

```RAPCOS.py``` is used to perform RAPCOS with loss function measuring the distance between natural vector and target natural vector. 

### Examples for RAP with Convex Hull

```
python RAP_Convex.py
```

where the ```main()``` function can be performed directly. 

```RAP_Convex.py``` is used to perform RAP with loss function measuring the distance from natural vector to convex hull. 

### Matlab Code for PCA program to obtain vertices

```NV.m```, ```pca_program.m``` and ```main.m``` are designed to perform principal component analysis to get vertices of hyperpyramid. (Version: MATLAB R2016a)

```inhull.m``` is a default function in Matlab, which is designed to check whether a point is inside a convex hull.

## Notice

The EPOCH parameter is the number of iteration to perform RAP and RAPCOS. The parameter EPOCH should not be set too large because when the loss approaches 0, it will be harder to converge. 

If the algorithm is stuck at some point, it needs to be stopped and run again, where some local minima might be achieved. 

### Authors

[Ruzhang Zhao](http://ruzhangzhao.com), Department of Biostatistics, Bloomberg School of Public Health, Johns Hopkins University

Shaojun Pei, Department of Mathematical Science, Tsinghua University

[Stephen Shing-Toung Yau](http://homepages.math.uic.edu/~yau/), Department of Mathematical Science, Tsinghua University

License

This project is licensed under the MIT License - see the LICENSE.md file for details
# NewGenomeSeqDetection
# NewGenomeSeqDetection

# RDB
This is an R package implementing RDB test proposed in Wang (2021)

To install this package in R, run the following commands:

```R
library(devtools) 
install_github("lakerwsl/RDB")
```

Example usage:

```R
library(WUNT)

m=50
d=100 
P=matrix(runif(m*d),nrow=m,ncol=d)
Z=rep(0,m)
Z[1:(m/2)]=1
rdb(P,Z)
```

#### References
Shulei Wang.
<b>Robust Differential Abundance Test in Compositional Data.</b>
2021.
[<a href="https://arxiv.org/pdf/2101.08765.pdf">arxiv</a>]
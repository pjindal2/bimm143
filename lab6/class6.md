Class6 R functions
================
Priya Jindal
10/17/2019

# This is H1

This is my work from class6 in **BIMM143**

## This is heading size 2

### This is heading 3 etc.

``` r
plot(1:10)
```

![](class6_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

## Practice reading files (again…)

``` r
read.table("test1.txt", sep=",", header=TRUE)
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
read.table("test2.txt", sep="$", header=TRUE)
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
read.table("test3.txt")
```

    ##   V1 V2 V3
    ## 1  1  6  a
    ## 2  2  7  b
    ## 3  3  8  c
    ## 4  4  9  d
    ## 5  5 10  e

Insert a code block: ctrl + alt + i

``` r
add <- function(x, y=1){
  #sum input x and y
  x + y
}
```

``` r
add(1)
```

    ## [1] 2

``` r
add(5,5)
```

    ## [1] 10

``` r
add(c(1,2,3),4)
```

    ## [1] 5 6 7

``` r
add(x=1, y=2)
```

    ## [1] 3

``` r
rescale <- function(x){
  rng <- range(x)
  (x-rng[1])/(rng[2] -rng[1])
}
```

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

# Section 2: of hands-on sheet

Install the **bio3d** package for sequence and structure analysis

installed once already

``` r
library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")

s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b

plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](class6_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](class6_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](class6_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

# Answer to questions

1.  read.pdb Returns a list of class “pdb”

2.  
<!-- end list -->

``` r
analyzeProtein <- function(pName){
  s1 <- read.pdb(pName)
  s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
  s1.b <- s1.chainA$atom$b
  plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
}

analyzeProtein("4AKE")
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): C:
    ## \Users\priya\AppData\Local\Temp\RtmpstpdoQ/4AKE.pdb exists. Skipping
    ## download

![](class6_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Homework Lecture 6 Problem 6
================
Priya Jindal
10/23/2019

``` r
library(bio3d)

# The input to this function is the kinase name that we want to study
# This function analyzes the protein and drug interation, to use it
# just call the function and pass it the kinase name as an argument
# The output of this function will be a graph showing the B-factor 
# trends for each kinase
analyzeProtein <- function(pName){
  s1 <- read.pdb(pName)
  s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
  s1.b <- s1.chainA$atom$b
  plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
}

#method call for kinase with drug
analyzeProtein("4AKE")
```

    ##   Note: Accessing on-line PDB file

![](PriyaJindalHw6Problem6_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
#method call for kinase with no drug
analyzeProtein("1AKE")
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

![](PriyaJindalHw6Problem6_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
#method call for kinase with drug
analyzeProtein("1E4Y")
```

    ##   Note: Accessing on-line PDB file

![](PriyaJindalHw6Problem6_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

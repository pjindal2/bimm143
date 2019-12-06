Structural Bioinformatics
================

``` r
# Read CSV 

data <- read.csv("Data Export Summary.csv")
data
```

    ##   Experimental.Method Proteins Nucleic.Acids Protein.NA.Complex Other
    ## 1               X-Ray   131278          2059               6759     8
    ## 2                 NMR    11235          1303                261     8
    ## 3 Electron Microscopy     2899            32                999     0
    ## 4               Other      280             4                  6    13
    ## 5        Multi Method      144             5                  2     1
    ##    Total
    ## 1 140104
    ## 2  12807
    ## 3   3930
    ## 4    303
    ## 5    152

Q1: Download a CSV file from the PDB site (accessible from “Analyze” -\>
“PDB Statistics” \> “by Experimental Method and Molecular Type”. Move
this CSV file into your RStudio project and determine the percentage of
structures solved by X-Ray and Electron Microscopy. Also can you
determine what proportion of structures are protein?

``` r
# Proportion of structures total
round((data$Total) / sum(data$Total) * 100, 2)
```

    ## [1] 89.07  8.14  2.50  0.19  0.10

``` r
#Proportion of protein
round(sum(data$Proteins) / sum(data$Total) * 100,2)
```

    ## [1] 92.71

``` r
# Here we are reading the 1HSG PDB structure and selecting teh protein component and 
# writing out only the protein component. Then we will repeat the same thing for 
# the ligand (drug) 

library(bio3d)
data <- read.pdb("1hsg.pdb")
protein <- atom.select(data, "protein", value = TRUE)
write.pdb(protein, file="1hsg.protein.pdb")
ligand <- atom.select(data, "ligand", value = TRUE)
write.pdb(ligand, file="1hsg.ligand.pdb")
```

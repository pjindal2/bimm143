Lab 12 Structural Bioinformatics -Drug Discovery
================
Priya Jindal
11/7/2019

\#\#Prep for docking

We want to produce a protein-only PDB file and a drug only PDB file.

\#\#Obtaining and inspecting our output structure

``` r
library(bio3d)

#dowload pdb file
file.name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

Read this PDB file into R so we can prepare it for further analysis.

``` r
hiv <- read.pdb("1hsg.pdb")

#GET THE PROTEIN COMPONENT OF THE MOLECULE
protein <- atom.select(hiv, "protein", value = TRUE)
write.pdb(protein, file="1hsg_protein.pdb")

#GET THE DRUG COMPONENET (LIGAND)
ligand <- atom.select(hiv, "ligand", value= TRUE)
write.pdb(ligand, file="1hsg_ligand.pdb")
```

## AutoDock Tools

Using AutoDock tools open the 1hsg\_protein.pdb file and use ADT to
assign charges and aton type to each atom in the protein.Next repeat
this with the ligand file.

Next set up a config files for the docking calculations based on the 3D
serach space where the ligand docking will be attempted (use the 3D
image of thr protein to determine this). Then use AutoDock Vina software
to dock the ligand and get the output of all the docking modes and a
table of affinities based on AutoDock Vina’s scoring function.

The best docked mode is the first entry in the all.pdbqt file.

\#Inspecting your docking results In order to visualize the docks and
compare to the crystal conformation of the ligand we will process the
all.pdbqt to a PDB format file that can be loaded into VMD. To do this
we will use R and the Bio3D package.

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

To assess the results quantitatively we will calculate the RMSD (root
mean square distance) between each of the docking results and the known
crystal structure using the bio3d package

Read the original ligand with added hydrogens that you produced earlier
and use the rmsd() function to compare to your docking results.

``` r
ori <- read.pdb("ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.649  4.206 11.110 10.529  4.840 10.932 10.993  3.655 10.996 11.222
    ## [11] 10.567 10.372 11.019 11.338  8.390  9.063  8.254  8.978

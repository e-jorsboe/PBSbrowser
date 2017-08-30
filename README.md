# PBSbrowser

The PBS browser uses population specific frequencies and size of each population,
for calculating Reynold's Fst both unweighted and weighted and PBS.

It also calculates gene Fst or PBS, which is a weigthed Fst calculated for a window spanned by the gene.


HOW TO USE:

1. Copy github to a new folder with your shiny-name
2. Create file called folderName with the name of the folder (shiny-name)
3. Run prepData.R, from its directory, with either .frq, bim and fam files from plink or .beagle, .fopt (+ .filter) and .qopt
files from NGSadmix.
4. Enjoy!
# PBSbrowser

The PBS browser uses population specific frequencies and size of each population,
for calculating Reynold's Fst both unweighted and weighted and PBS.

It also calculates gene Fst or PBS, which is a weigthed Fst calculated for a window spanned by the gene.

These R packages are reqiured for now: (02-10-2017)
Should install 'shiny', 'inline', 'data.table' and 'parallel' R packages

HOW TO USE:

1. Copy github to a new folder with your shiny-name
(1.5. Add password as 'passWord' variable in server.R script in folder)
2. Run prepData.R, from its directory, with either .frq, bim and fam files from plink or .beagle, .fopt (+ .filter) and .qopt from NGSadmix. Or .bim file, Q and P files from
ADMIXTURE. For more info run 'Rscript prepData.R'
3. Either run locally using runApp() from the shiny package, or run on server and access shiny via browser!
4. Enjoy!
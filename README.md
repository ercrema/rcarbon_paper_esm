# Supplementary Material for "Inference from large sets of radiocarbon dates: software and methods"
This repository contains all data and scripts required to fully reproduce all analyses presented in the following paper:

Crema, E.R., Bevan,A., 2020. [Inference from large sets of radiocarbon dates: software and methods](https://doi.org/10.1017/RDC.2020.95). Radiocarbon, DOI: 10.1017/RDC.2020.95

# File Structure
- `/data/*.csv` ... Contains CSV files of sample case studies
- `runFigures.R` ... R script for generating figures in the paper
- `*.tiff` ... Tiff version of the figures in the paper
- `esm_rcarbon_paper.Rmd` ... Rmarkdown containing notes and scripts for generating manuscript and supplementary figures.
- `esm_rcarbon_paper.pdf` ... Rendered version of `esm_rcarbon_paper.Rmd`


# R Settings
```
R version 4.0.2 (2020-06-22)

attached base packages:
[1] stats     graphics  grDevices utils     methods   base     

other attached packages:
[1] maptools_1.0-1 sp_1.4-2       rcarbon_1.4.0 

loaded via a namespace (and not attached):
 [1] lattice_0.20-41       codetools_0.2-16      deldir_0.1-28        
 [4] snow_0.4-3            foreach_1.5.0         grid_4.0.2           
 [7] nlme_3.1-148          tensor_1.5            spatstat_1.64-1      
[10] rpart_4.1-15          goftest_1.2-2         Matrix_1.2-18        
[13] startup_0.14.1        spatstat.utils_1.17-0 splines_4.0.2        
[16] iterators_1.0.12      tools_4.0.2           foreign_0.8-80       
[19] polyclip_1.10-0       spatstat.data_1.4-3   abind_1.4-5          
[22] xfun_0.16             parallel_4.0.2        compiler_4.0.2       
[25] mgcv_1.8-31           doSNOW_1.0.18         knitr_1.29

```
# Funding
This research was partly funded by the ERC grant Demography, Cultural Change, and the Diffusion of Rice and Millets during the Jomon-Yayoi transition in prehistoric Japan (ENCOUNTER) (Project N. 801953, PI: Enrico Crema).

# Licence
CC-BY 3.0

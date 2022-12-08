This git includes all processed data that have been generated and used in the study "Random genetic drift sets an upper limit on mRNA splicing accuracy in metazoans", as well as the scripts used to analyze the data and to generate the figures.

Developed by Florian Bénitière, Anamaria Necsulea and Laurent Duret.
Université de Lyon, Université Lyon 1, CNRS, Laboratoire de Biométrie et Biologie Évolutive UMR 5558, F-69622 Villeurbanne, France.


### Directories content ###

- The directories "figure generator" and "figure_supp generator" contain R scripts used to produce figures stored in the respective directories "figure" and "figure_supp".

- The directory "table_supp" contain Rmarkdown scripts used to generate supplementary tables, which are stored in the same directory

- The directory "data" contain processed data under table format used to generate figures and analyses

- Stable images used in pannels are stored in "images_library" and intermediate pannels generated by R scripts to produce figures are stored in the directory "pannels"

- The "pipelines" directory contain the pipeline to calculate the dN/dS and the one used to analyse various alternative splicing characteristics as part of our study


###
Some packages that you may encounter in these R scripts should be install:
- stringr
- kableExtra
- cowplot
- ggplot2
- imager
- RColorBrewer
...


###
The script "values_in_text_paper_generator.R" retrieves the values found in the main text of the paper from the data

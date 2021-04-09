# Genomic-epidemiology-of-SARS-CoV2-in-Russia
This repository contains scripts and data used in Genomic epidemiology of SARS-CoV-2 in "Russia reveals recurring cross-border transmission throughout 2020" paper
Here is a brief description of all files:

Folder "TransmissionLineages" contains files with russian singletons, stem clusters and transmission lineages as in (Komissarov et al, 2021) paper (https://www.nature.com/articles/s41467-020-20880-z). It also contains file with Pangolin clades for each russian sample (full_cat_pangolin_clades.out) as well as R script that uses all these files to make figure 3 from the "Genomic epidemiology of SARS-CoV-2 in Russia reveals recurring cross-border transmission throughout 2020" Manuscript.

Folder "Map" contains table with region, collecting date, Pangolin lineage and russian transmission lineage for each russian sample. It alsop contains an R script draw_maps.PangilinColors.Figures1_and_4.R that uses this table to create Figures 1 and 4 from the Manuscript.

Folder "Cross-bordeerTransmissions" contains Perl scripts that find cross-border transmissions inside (russian_lineages_from_nonrussian_dates.strict.with_seq_number.with_countries.pl) and out of (nonrussian_lineages_from_russian_dates.strict.with_seq_number.with_countries.pl) Russia using phylogenetic trees for main SARS-CoV2 clades (files .nwk). It also contains results of node dating with node.dating procedure from ape R package (files .nodes_dates_old_labels.debugged.csv) for each tree. File Avia_and_Sequences.txt is a file with information about air passenger traffic and SARS-CoV2 sampling intensity in Russia for each month. Script Timeline_ImportsExports.Dates.Figure5.R makes Figure 5 from the manuscriopt using these files. 

# Orthogroups

A little vagueley defined right now, but code related to processing the results from orthofinder

## To do
* read the list of single copy orthologs
* read the list of unassigned genes

how shall we store these? sc orthologs can probably be an added field in the orthogroup class
unassigned genes doesn't fit well, maybe need another class to store the genes in each proteome with assigned orthogroup

* how to store species metadata?
* select species with large differences in number or proteins based on metadata contrasts

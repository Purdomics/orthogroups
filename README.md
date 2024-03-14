# Orthogroups
## outliers.py
Compare a target and non-target group and identigy orthogroups that are expanded or contracted in the target group
relative to the not target group. Sequence counts in each OG are converted to a standard normal deviate using the
50% trimmed mean and standard deviation
```
usage: outliers.py [-h] [-g ORTHOGROUP] [-n NTOP] [-f FRACTION] [-t TARGET] [-j JSON] [-v TSV]
Find expanded and contracted orthogroups
optional arguments:
  -h, --help                            show this help message and exit
  -g ORTHOGROUP, --orthogroup ORTHOGROUP
                                        Orthogroups.tsv file
  -n NTOP, --ntop NTOP                  number of top/bottom groups to report
  -f FRACTION, --fraction FRACTION      total fraction of observations to trim in mean/std.dev. calculation
  -t TARGET, --target TARGET            comma delimited string with list of target organisms
  -j JSON, --json JSON                  JSON output file
  -v TSV, --tsv TSV                     file for normalized counts in TSV format
```
## og_interpro.py
Submit each sequencea in a list of OGs to interproscan. Results are saved as pickled json.
```
usage: og_interpro.py [-h] [-g ORTHOGROUP] [-o OUT] [-s]

Run interpro on selected orthogroups

optional arguments:
  -h, --help                            show this help message and exit
  -g ORTHOGROUP, --orthogroup ORTHOGROUP
                                        Orthogroups file
  -o OUT, --out OUT                     output directory for results
  -s, --skip                            skip (default: False)
```

## og_analyze_ips_result.py
```
usage: og_analyze_ips_result.py [-h] [inputfilename] [outputfilename]

Analyze interpro result for orthogroups

positional arguments:
  inputfilename   Intput file name (default: pklfiles/*.pkl)
  outputfilename  Output file name (default: og.analysis.txt)

optional arguments:
  -h, --help      show this help message and exit
```
* example output
  
```
og_analyze_ips_result.py 2024-03-14 14:48:55
        OG input: pklfiles/*.pkl
        Output directory: og.analysis.txt

OG0000525       99 sequences    
        json unreadable or no matches pklfiles\OG0000525_jgi_Cap6580_1_796925_estExt_Genewise1Plus.C_2840032.pkl
        json unreadable or no matches pklfiles\OG0000525_jgi_Cap6580_1_842022_fgenesh1_kg.61_#_80_#_TRINITY_DN7932_c0_g1_i1.pkl

        Feature: SM00829: PKS_ER_names_mod (IPR: IPR020843 - Polyketide synthase, enoylreductase domain)
        Number: 80      High score: 0.0078      Mean score: 7e-05       Low score: 3.1e-07
        Number: 80      GO:0016491: (MF) oxidoreductase activity

        Feature: PF00107: ADH_zinc_N - Zinc-binding dehydrogenase (IPR: IPR013149 - Alcohol dehydrogenase-like, C-terminal)
        Number: 91      High score: 2.5e-06     Mean score: 3.4e-16     Low score: 1.3e-22

        Feature: PF08240: ADH_N - Alcohol dehydrogenase GroES-like domain (IPR: IPR013154 - Alcohol dehydrogenase-like, N-terminal)
        Number: 91      High score: 0.0041      Mean score: 9.1e-16     Low score: 6.2e-21

        Feature: cd08297: CAD3 - CAD3 (IPR: None - )
        Number: 15      High score: 2.294e-114  Mean score: 3.4e-120    Low score: 4.53331e-128

        Feature: PF00400: WD40 - WD domain, G-beta repeat (IPR: IPR001680 - WD40 repeat)
        Number: 1       High score: 740.0       Mean score: 7.4e+02     Low score: 740.0
        Number: 1       GO:0005515: (MF) protein binding

        Feature: SM00320: WD40_4 (IPR: IPR001680 - WD40 repeat)
        Number: 1       High score: 0.044       Mean score: 0.044       Low score: 0.044
        Number: 1       GO:0005515: (MF) protein binding

        Feature: PF08297: U3_snoRNA_assoc - U3 snoRNA associated (IPR: IPR013268 - U3 snoRNA associated)
        Number: 1       High score: 1.9e-09     Mean score: 1.9e-09     Low score: 1.9e-09
        Number: 1       GO:0006364: (BP) rRNA processing
        Number: 1       GO:0030515: (MF) snoRNA binding

        Feature: PS00059: ADH_ZINC - Zinc-containing alcohol dehydrogenases signature. (IPR: IPR002328 - Alcohol dehydrogenase, zinc-type, conserved site)
        Number: 1       High score: 1   Mean score: 1   Low score: 1
        Number: 1       GO:0016491: (MF) oxidoreductase activity
        Number: 1       GO:0008270: (MF) zinc ion binding
```
## To do
* add metadata to allow more flexible selection by og_outliers.py
  * how to store species metadata?

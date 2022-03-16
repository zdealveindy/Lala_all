[![DOI](https://zenodo.org/badge/464408388.svg)](https://zenodo.org/badge/latestdoi/464408388)

# Lalashan elevation transect - all data

## Paper
This repository holds the code and data accompanying the following paper:

"Climate and soil differentially affect species, trait and diversity patterns of woody overstory and fern understory in a subtropical forest along an elevation gradient in Taiwan" [2022], by Kenny Helsen, Yen-Cheng Shen, Tsung-Yi Lin, Chien-Fan Chen, Chu-Mei Huang, Ching-Feng Li & David Zelený

> *Abstract*:  Questions. Although the relative importance of climate in abiotic filtering is higher for woody than herbaceous species assemblages, it is unclear whether this pattern is also reflected between the woody overstory and herbaceous understory of forests. The understory might respond more to small-scale soil variation, next to experiencing additional abiotic filtering through overstory effects on light and litter quality. We explored the proportional importance of climate and soil on the species, trait and (functional) diversity patterns of both the forest overstory and fern and lycophyte understory.
Location. Subtropical forest along an elevational gradient from 850 to 2100 m a.s.l. in Northern Taiwan.
Methods. We measured nine functional traits expected to respond to soil nutrient or climatic stress for woody overstory species and understory ferns and lycophytes. Next, we performed parallel constrained ordinations on over- and understory species and trait composition, and multiple regression for species and functional diversity, using measured climate proxies and soil variables as predictors.
Results. Climate was more important than soil in predicting the species composition of both vegetation layers and trait composition of the understory. The stronger than expected effect of climate for the understory was likely due to fern and lycophytes' higher vulnerability to drought, while the higher importance of soil for the overstory trait composition seemed driven by deciduous species. The environmental drivers affected different response traits in both vegetation layers, and the overstory had additional effects on understory traits, resuling in a disconnection of community-level trait values across layers. Interestingly, species and functional diversity patterns could be almost exclusively explained by climate effects for both layers. 
Conclusions. This study illustrates that abiotic filtering can differentially affect species, trait and diversity patterns and can be highly divergent for forest overstory and fern understory vegetation, and should consequently not be extrapolated across vegetation layers.

## Code
The file `Lala_woody_ferns_spe_traits_diversity.R` contains complete R code (with comments) used to conduct all analyses and draw all figures.

## Data
All data are provided as tab-delimited *.txt files:
- `lala_100_spe.txt` - species composition matrix (including both woody and fern species);
- `lala_100_env.txt` - matrix with environmental variables;
- `lala_woody_traits.txt` - measured leaf traits of all woody species in the analysis;
- `lala_fern_traits.txt` - measured leaf traits of all fern and lycophyte species in the analysis;
- `lala_100_checklist.txt` - checklist of species, including abbreviations.

## Persistent archiving
This repository is mirrored on [Zenodo](https://zenodo.org/), providing long-term access. Please consider citing the accompanying paper if you re-use this code for academic purposes.

## License
[![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg



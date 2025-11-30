
# Haruka


**A Python package for identifying salient and background spatial domains in large-scale spatial omics analysis.**

Haruka identifies salient and background spatial domains (DSEP) across multi-slice datasets and provides downstream analyses based on the discovered domains.

</p>
<p align="center">
  <img src="https://github.com/C0nc/Haruka/blob/main/main_figs/image.png" width="800px">
</p>

## Getting started


Please refer to the  
- MERIFSH simulation [Tutorial][link-tutorial_1] (Can be downloaded by ) 
- MERISH brain aging dataset [Tutorial][link-tutorial_2] (Can be downloaded via [link][aging_data])

## Installation

1. Create a conda environment
```bash
conda env create -f environment.yml --name haruka
conda activate haruka
```
2. Install the modified scVI-tools dependency
```bash
cd scvi-tools
pip install -e ".[dev]"
```

3. Install the mini-batch version scENVI COVET calculation
```bash
pip install 'scenvi==0.4.5' --no-deps
```

(Optional) 

Install MaxFuse [link][maxfuse] and scSLAT [link][scSLAT] for reproduce the downstream analysis


## Acknowledgement:

This project is built based on the codebase of scVI [link][scVI]. cellcharter [link][cellcharter] and contrastiveVI [link][contrastiveVI]. 

## Contribution

If you found a bug or you want to propose a new feature, please use the [issue tracker][issue-tracker].

[issue-tracker]: https://github.com/C0nc/Haruka/issues
[link-tutorial_1]: https://github.com/C0nc/Haruka/blob/main/repo_sim.ipynb
[link-tutorial_2]: https://github.com/C0nc/Haruka/blob/main/aging_reproduce.ipynb
[maxfuse]: https://github.com/shuxiaoc/maxfuse
[scSLAT]: https://github.com/gao-lab/SLAT.git
[cellcharter]: https://github.com/CSOgroup/cellcharter
[scVI]: https://scvi-tools.org/
[contrastiveVI]: https://github.com/suinleelab/contrastiveVI
[aging_data]: https://doi.org/10.5281/zenodo.13883177

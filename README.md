# confWorks

Convenient python-wrapper for handling multconformer RDKit molecules and conducting geometry optimization and conformational sampling routines (GFN-xTB semi-empirical theory levels).

Coded by Galymzhan Moldagulov (Grad student @ IBS-CARS & UNIST, KR)

You can easily:
- Handle RDKit Mols with multiple conformers: saving, reading files, storing and loading conformer specific attributes from SDfs;
- Calculating xTB energies and vibrational frequencies;
- Conduct geometry optimization (xTB) for all or only selected conformers of a Mol;
- Conduct conformer sampling (CREST) for all or only selected conformers of a Mol;
- Calculate Boltzmann weights for individual conformers;
- Visualize the conformational landscape of RDKit Mols;
- Inspect individual geometry optimization trajectories.

## Requirements
- Python (version ≥ 3.9);
- RDKit (version ≥ 2022);

## Installation

Install *xTB* separately via conda:
```conda install xtb``

Then install the *confWorks* package:
1. ```git clone https://github.com/moldagulovg/confWorks.git```
2. cd into cloned *confWorks* folder;
3. ```pip install -e .```
4. now you may import to your python code ```import confworks```



## Short Overview (TBA):
- read_multiconf_sdf;
- save_multiconf_sdf;
- xtb_SP;
- optimize_molecule;
- conformer_search;

![](https://moldagulovg.github.io/confWorks/docs/assets/quinine_2D.jpg)
![](https://moldagulovg.github.io/confWorks/docs/assets/quinine_3D.jpg)
![](https://moldagulovg.github.io/confWorks/docs/assets/quinine_3D_ensemble.jpg)
![](https://moldagulovg.github.io/confWorks/docs/assets/rmsd_matrix.jpg)
![](https://moldagulovg.github.io/confWorks/docs/assets/rmsd_hist.jpg)
![](https://moldagulovg.github.io/confWorks/docs/assets/conf_space_tsne.jpg)
![](https://moldagulovg.github.io/confWorks/docs/assets/conf_space_energy_landscape.jpg)

## What's next?
- Further updates will include conformer sampling using CREST;
- xTB metadynamics simulations.


## References:
1. RDKit: Open-source cheminformatics. (https://www.rdkit.org, accessed on 21.7.2025), https://doi.org/10.5281/zenodo.7671152;
2. C. Bannwarth, E. Caldeweyher, S. Ehlert, *et al.* Extended tight-binding quantum chemistry methods. *WIREs Comput. Mol. Sci.* **2021**; 11:e1493, https://doi.org/10.1002/wcms.1493;
3. C. Bannwarth, S. Ehlert, S. Grimme, *J. Chem. Theory Comput.*, **2019**, 15, *3*, 1652–1671, https://pubs.acs.org/doi/10.1021/acs.jctc.8b01176;
4. P. Pracht, F. Bohle, S. Grimme, *Phys. Chem. Chem. Phys.*, **2020**, 22, *14*, 7169-7192, https://doi.org/10.1039/C9CP06869D;
5. S. Spicher, S. Grimme, *Angew. Chem. Int. Ed.*, **2020**, 59, 15665, https://doi.org/10.1002/anie.202004239;
6. L. Turcani, A. Tarzia, F.T. Szczypiński, K.E. Jelfs, *J. Chem. Phys.*, **2021**, 154 (21): 214102, https://doi.org/10.1063/5.0049708;
7. G. Fraux, R.K. Cersonsky, M. Ceriotti, *Journal of Open Source Software*, **2020**, 5, *51*, 2117, https://doi.org/10.21105/joss.02117.


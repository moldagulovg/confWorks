# confWorks
Convenient workflow for generation of geometry optimized conformers on GFN-xTB semi-empirical theory levels.

Coded by Galymzhan Moldagulov (Grad student @ IBS-CARS & UNIST, KR)

You can easily:
- explore chemical space in the interactive Faerun map;
- search complexes by their CSD code;
- switch between color schemes in the plot legend;
- look at 3D structures of individual complexes 

![](https://moldagulovg.github.io/rdkit-xtb-geomopt/docs/assets/quinine_2D.jpg)
![](https://moldagulovg.github.io/rdkit-xtb-geomopt/docs/assets/quinine_3D.jpg)
![](https://moldagulovg.github.io/rdkit-xtb-geomopt/docs/assets/quinine_3D_ensemble.jpg)
![](https://moldagulovg.github.io/rdkit-xtb-geomopt/docs/assets/rmsd_matrix.jpg)
![](https://moldagulovg.github.io/rdkit-xtb-geomopt/docs/assets/rmsd_hist.jpg)
![](https://moldagulovg.github.io/rdkit-xtb-geomopt/docs/assets/conf_space_tsne.jpg)
![](https://moldagulovg.github.io/rdkit-xtb-geomopt/docs/assets/conf_space_energy_landscape.jpg)

## Methods:
- Sampled 10-108k mononuclear organometallic complexes from tmQM dataset and CSD <sup>[1, 2]</sup>;
- Vectorized molecules were embedded/projected into 2D TMAP representation <sup>[3]</sup>;
- Molecules were vectorized by MHFP (MinHash Fingerprint)<sup>[4]</sup>;
- Chemical space was visualized with interactive Faerun module <sup>[5]</sup>;
- Visualization of individual molecules in various 3D representations was impleted using 3Dmol.js <sup>[6]</sup>.

## Further updates will include conformer sampling using CREST, as well as in line use of ORCA for QM calculations.

## References:
1. *J. Chem. Inf. Model*. **2020**, 60, *12*, 6135–6146;
2. *Acta Cryst.* **2016**, B72, 171-179;
3. *J. Cheminform.*, **2018**, 12, *12*;
4. *J. Cheminform.*, **2018**, 10, *66*;
5. *Bioinformatics*, **2018**, 34, *8*, 1433–1435;
6. *Bioinformatics*, **2015**, 31, *8*, 1322–1324.

# LBM-PSM-DEM-Extension
This project extends the LBM-DEM coupling interface by Palabos, LIGGGHTS, and Philippe Seil, addressing cases where a single LBM lattice intersects multiple particles and resolving issues with particles exiting the computational domain. It integrates the DLVO collision model, enabling realistic particle interactions at micro- and nano-scales.
For specific installation procedures, please refer to https://github.com/ParticulateFlow/LBDEMcoupling-public.
Simply copy the code from src-PSM in this repository into LBDEM-coupling/src, and copy the code from src-DEM into LIGGGHTS/src.
It should be noted that after modifying LIGGGHTS/src, the constitutive codes must be added to the Makefile.lib file, and liggghts must be recompiled.

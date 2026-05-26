# LBM-PSM-DEM-Extension
## About
This project extends the LBM-DEM coupling interface by Palabos, LIGGGHTS, and Philippe Seil, addressing cases where a single LBM lattice intersects multiple particles and resolving issues with particles exiting the computational domain. It integrates the DLVO collision model, enabling realistic particle interactions at micro- and nano-scales.
## Installation
For specific installation procedures, please refer to https://github.com/ParticulateFlow/LBDEMcoupling-public.
Simply copy the code from src-PSM in this repository into LBDEM-coupling/src, and copy the code from src-DEM into LIGGGHTS/src.
It should be noted that after modifying LIGGGHTS/src, the constitutive codes must be added to the Makefile.lib file, and liggghts must be recompiled.
## Reference
Wang, S., Hou, P., Liu, Q., Sang, G., Liang, X., Dou, F., et al. (2025). Microparticle transport and clogging mechanisms in 3D complex rock fractures based on coupled lattice Boltzmann and discrete element method simulations. International Journal of Rock Mechanics and Mining Sciences, 195, 106259. https://doi.org/10.1016/j.ijrmms.2025.106259
Wang, S., Hou, P., Wu, Z., Liu, Q., Sang, G., Rabatuly, M., et al. (2026). Mechanisms of suspended microparticle clogging in
constricted channels: Insights from a lattice Boltzmann‐discrete element simulation. Water Resources Research, 62, e2025WR042425. https://doi.org/10.1029/2025WR042425

# MDistEn
Algorithm for the computation of Multichannel Distribution Entropy (MDistEn)

This code implements the MDistEn algorithm for the computation of the Distribution Entropy value in multichannel systems, as described in:

[1] A. Gargano, M. Nardelli, E.P. Scilingo, "Exploring Multivariate Dynamics of Emotions Through Time-Varying Self-Assessed Arousal and Valence Rating", 

MDistEn relies on a novel method for the reconstruction of multivariate phase spaces, in which each time series is embedded using its proper time delay, following a similar approach used for the MCI index (https://doi.org/10.1016/j.physa.2019.121543). However, MDistEn uses the information from the phase space distances to measure the spatial complexity of the multivariate trajectory and investigates the complex dynamics of the trajectory.

The algorithm was used on continuously annotated time series, as described in the article.

_________________________________________________________________________

Copyright (C) 2024 Andrea Gargano, Mimma Nardelli

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

If you use this program in support of published research, please include a citation of the reference [1]. If you use this code in a software package, please explicitly inform the end users of this copyright notice and ask them to cite the reference above in their published research.
__________________________________________________________________________

After downloading all available functions, to use this software in Matlab you have to call the MDistEn function in your working path/folder. Type 'help MDistEn' from the command window to receive help on the command's syntax and input/output arguments.

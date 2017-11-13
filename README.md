# BTW_SandpileSimulation
Implementation of the sandpile model from [Bak–Tang–Wiesenfeld, 1987]. 

*As default the code drops 15 grains on a random position as BTW,1987.  To evidence pattern formation, 1 grain is dropped at the center of the plane for N times. 

BTW model:
Self-organized criticality: An explanation of the 1/f noise. https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.59.381
More info: http://guava.physics.uiuc.edu/~nigel/courses/563/Essays_2012/PDF/banerjee.pdf

Tags: Self emerging patterns, Power Law.

Files:

sandpile_jb.py:    Plot an png image of the sandpile in each time step as in BTW,1987.
sandpile_jb_pattern.py:    Plot an png image of the sandpile in each time step to evidence the pattern formation.  One grain is dropped every step in the center of the plane.
sandpile_jb_v2.py: Export a csv file containing data of the number of grains dropped in each avalanche to further build an histogram. 

Results:
![alt text](https://github.com/JoseBarreiros/BTW_SandpileSimulation/blob/master/Media/avalanches.png)

Pattern formation:
![alt text](https://github.com/JoseBarreiros/BTW_SandpileSimulation/blob/master/Media/avalanches.gif)

Avalanches count vs Lifetime:
![alt text](https://github.com/JoseBarreiros/BTW_SandpileSimulation/blob/master/Media/histo.png)



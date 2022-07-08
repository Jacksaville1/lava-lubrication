# Lava Flow Modelling
This repository contains the code used to generate the results in the paper "Predicting Safe Zones within Lava Flows over Topography" by Saville, Hinton and Huppert.
The paper uses lubrication theory (a shallow-flow approximation) to model isoviscous lava flows over topography. Time-dependent lava simulations are produced by solving the governing equation

$ \frac{\partial h}{\partial t} = \frac{\rho g}{3 \mu} \nabla \cdot [h^3 \nabla (E+h)] $,    

while steady state simulations are produced by setting $\frac{\partial h}{\partial t} =0 $. This gives the steady governing equation 

$ \nabla \cdot [h^3 \nabla (E+h)]  =0 $.

In section 3 of the paper, this steady governing equation is solved for flow over mathematically idealised topographic features on inclined planes using 'steadysolver.m'. This program requires a topographic feature be specified, as well as the background slope, $\beta$, and lava thickness at the upstream line source, $h_\infty$. The program uses an adaptive solver from MATLAB's Partial Differential Equation Toolbox, and the upstream depth is iteratively reduced to $h_\infty$. The final iteration of the solver increases the mesh resolution to produce a high-resolution figure.

In section 4 of the paper, we turn our attention to real-world topography and model the transient flow at Marcath Volcano, Nevada in the eruption 35kyr ago. The lava field produced in this eruption is well preserved, and has been measured using LiDAR by Younger et al. (2019). The outline of the lava field from Younger et al. (2019) is contained in the file 'lavafield.tif', while the DEM, obtained from SRTM is contained in the file 'topography.tif'. These files are read by the program 'import_marcath_data.m', which constructs an approximation of the pre-eruptive palaeotopography at Marcath by subtracting the average lava thickness (11m) from the DEM. This is smoothed using a Gaussian filter with radius $R$
, which is specified. 
The transient solution is obtained using the program 'transientsolver.m'. This takes inputs of topography, eruption duration, effusion rate, 
$q$, 
lava viscosity, 
$\mu$ 
, and lava density, 
$\rho$
. The program calculates lava thickness at every mesh point for all times between the eruption start and end. These data are stored in variables 
$hsol$ 
and 
$t$.




References:

Younger, Z.P., Valentine, G.A. & Gregg, T.K.P. (2019) A’ā lava emplacement and the significance of rafted pyroclastic material: Marcath volcano (Nevada, USA). Bulletin of Volcanology, 81(9), doi:10.1007/s00445-019-1309-6.

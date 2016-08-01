# electro-optics

A variety of Matlab functions to calculate electrorefraction and electroabsorption effects in semiconductors.

These functions are for calculating the photon energy-dependent shifts in the real (n) and imaginary (k) parts of the refractive index of semiconductors, due to carrier induced effects (bandfilling, bandgap shrinkage and free carrier absorption) and electric field effects (the Franz-Keldysh effect). Some of them return the change in absorption coefficient (alpha, in m^-1) and some the change in the imaginary part of the refractive index (k). They are based on the models of Bennett [1,2] and Weber [3]. The function KKtrans.m uses the Maclaurin method to perform the numerical itegration, as described in Ohta [4]. Please refer to these publications to fully understand the applicability and limitations of these models.

The script FK_and_carrier_Dn_Dk.m calculates the total shift in n and k due to all of these effects, in a material with a set of electron concentrations N, hole concentrations P and electric fields Q, at a range of photon energies E (see the comments in the script for details of other material parameters required). Note, to perform the Kramers Kronig transform, the change in absorption is required over a large range of photon energies. Therefore please use a reasonably large range for E, even if only a limited number of photon energies are required (this is something that could be addressed in subsequent versions). 

The Franz Keldysh effect and free carrier absorption sections of the script can be run on their own, provided that all the required parameters are in the workspace. However, the bandgap shrinkage and bandfilling sections should always be run one after the other as the bandfilling calculation uses the modified bandgap, taking bandgap shrinkage into account. 

The script returns two arrays, DN and DK, of dimensions (length(N),length(E)), where the rows correspond to different carrier concentrations and the columns to different photon energies, for example DN(1,1:end) is a row vector containing the photon energy-dependent total shifts in refractive index for carrier concentrations N(1) and P(1) and electric field Q(1). Each section also returns an array of the same dimensions, containing the shifts in n and k due to the individual effect being calculated in that section.

Future improvements I would like to make:
- addition of the quantum confined Stark effect 
- addition of the linear and quadratic electro-optic effects
- modify the script so that it's possible to read in successive sets of material parameters from a database (e.g. so that shifts for an entire layer structure can be calculated in one go)
- check the results against up-to-date experimental data




References

[1] B. R. Bennett, R. A. Soref, and J. A. Del Alamo, “Carrier-induced change in refractive index of InP, 
GaAs and InGaAsP,” IEEE Journal of Quantum Electronics, vol. 26, no. 1, pp.
113–122, Jan. 1990
[2] J.-P. Weber, “Optimization of the carrier-induced effective index change in InGaAsP waveguides-
application to tunable Bragg filters,” IEEE Journal of Quantum Electronics, vol. 30, no. 8, 
pp. 1801–1816, Aug. 1994.
[3] B. R. Bennett and R. A. Soref, “Electrorefraction and electroabsorption in InP, GaAs, GaSb, InAs, 
and InSb,” IEEE Journal of Quantum Electronics, vol. 23, no. 12, pp. 2159–2166, Dec. 1987.
[4] K. Ohta and H. Ishida, “Comparison among Several Numerical Integration Methods for Kramers-Kronig Transformation,” Appl Spectrosc, vol. 42, no. 6, pp. 952–957, Aug. 1988.

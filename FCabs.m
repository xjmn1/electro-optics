function [DA,DN] = FCabs(N,P,E,ri,me)
%Function to calculate shift in absorption coefficient (DA) and refractive index (DN)
%due to free carrier absorption in semiconductors
%Based on equations 14-17 in J.-P. Weber, 
%“Optimization of the carrier-induced effective index change in InGaAsP waveguides-
%application to tunable Bragg filters,” 
%IEEE Journal of Quantum Electronics, vol. 30, no. 8, pp. 1801–1816, Aug. 1994.

%valid for wavelengths between 1 and 2 um, based on n-GaAs

%Input arguments:
%N = electron concentration (m^-3)(scalar)
%P = hole concentration (m^-3) (scalar)
%E = scalar or row vector of photon energies (eV)
%me = electron effective mass (m0)
%ri = refractive index of material: can be scalar or row vector with same
%size as E (if E is scalar, ri must be scalar)

hbar=1.055e-34; % reduced Planck's constant (kg m^2 s^-1)
e=1.602e-19; % electron charge (C)
eps_0=8.854e-12; % Vacuum permittivity (F m^-1)
c=2.998e8; % speed of light (m s^-1)

bE=E.*3.657;

%Electrons
DA_E=3.1e-22.*N.*ones(1,length(E)); 
DN_E=-(hbar.^2.*e.^2./(2.*eps_0.*ri.*me)).*(N./E.^2); 

%Holes
DA_H=4.252e-20.*exp(-3.657.*E).*P; 
DN_H=-(hbar.*c.*4.252e-20./(pi.*2.*E.*e)).*(exp(-bE).*(-expint(-bE)-1i.*pi)+exp(bE).*expint(bE)).*P; 

DA=DA_E+DA_H; %Total change in absorption coefficient (m^-1)
DN=DN_E+DN_H; %Total change in refractive index

end


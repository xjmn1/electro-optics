function [BGS] = BGshrink(me,eps_s,kappa,N,P,type)
% Function to calculate carrier-induced bandgap shrinkage (BGS) in semiconductors, using equations 17 and 18 in
% B. R. Bennett, R. A. Soref, and J. A. Del Alamo, 
% “Carrier-induced change in refractive index of InP, GaAs and InGaAsP,” 
% IEEE Journal of Quantum Electronics, vol. 26, no. 1, pp. 113–122, Jan.
% 1990

% BGS is the calculated bandgap shrinkage, in the same units as kappa

% Input arguments:
% N and P are scalars or row vectors of carrier concentrations in m^-3
% N and P must be same length if vectors
% me is electron effective mass  (m0)
% eps_s is the static dielectic constant of the material (eps_0)
% kappa is a fitting parameter (same units as bandgap energy - J or eV)
% type=an integer with value 1 or 2: for n-doping and electron-hole plasma injection, set type=1; for p-doping
% set type=2

   
Chi_cr=1.6e30.*(me./(1.4.*eps_s)).^3; 

if type==1
    
    BGS=-(kappa./eps_s).*((N./Chi_cr)-1).^(1/3); %bandgap shrinkage caused by excess electrons (J or eV) 
    BGS(N<Chi_cr)=0; %bandgap shrinkage = 0 for carrier concentration less than Chi_cr
    
elseif type==2
    
    BGS=-(kappa./eps_s).*((P./Chi_cr)-1).^(1/3); %bandgap shrinkage caused by excess holes (J or eV) 
    BGS(P<Chi_cr)=0;
    
end

end


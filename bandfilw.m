function [BF] = bandfilw(E,N,P,T,Eg,me,mhh,mlh,C)
%Function to calculate change in absorption coefficient (BF) due to bandfilling
%in semiconductors, using the model provided in B. R. Bennett, R. A. Soref, and J. A. Del Alamo, 
%“Carrier-induced change in refractive index of InP, GaAs and InGaAsP,” IEEE Journal of Quantum Electronics, 
%vol. 26, no. 1, pp. 113–122, Jan. 1990
%and modified in J.-P. Weber, “Optimization of the carrier-induced effective index change 
%in InGaAsP waveguides-application to tunable Bragg filters,” IEEE Journal of Quantum Electronics, 
%vol. 30, no. 8, pp. 1801–1816, Aug. 1994.

%Input arguments
%E = row vector of photon energies (J)
%N = electron concentration (m^-3) (scalar)
%P = hole concentration (m^-3) (scalar)
%T = temperature (K) (scalar)
%Eg =material bandgap (J) (scalar)
%me, mhh, mlh are electron, heavy hole and light hole effective masses (kg)
%(scalars)
%C is fitted parameter from absorption data (m^-1 J^1/2) (scalar) (see
%Bennett 1990)

kB=1.381e-23; % Boltzmann's constant (m^2 kg s^-2 K^-1)
hbar=1.055e-34; % reduced Planck's constant (kg m^2 s^-1)

mu_ehh=me*mhh/(me+mhh); % electron-heavy hole reduced effective mass (kg)
mu_elh=me*mlh/(me+mlh); % electron-light hole reduced effective mass (kg)
mdh=(mhh^(3/2)+mlh^(3/2))^(2/3);% density of states effective mass for holes (kg)

Chh=C*(mu_ehh^(3/2)/(mu_ehh^(3/2)+mu_elh^(3/2))); %contribution to C from heavy holes 
Clh=C*(mu_elh^(3/2)/(mu_ehh^(3/2)+mu_elh^(3/2))); %contribution to C from light holes

Nc=2*(me*kB*T/(2*pi*hbar^2))^(3/2); % effective density of states in CB (eqn 9a in Bennett 1990)
Nv=2*(mdh*kB*T/(2*pi*hbar^2))^(3/2); % effective density of states in VB (eqn 9b Bennett 1990)


Efc=kB.*T.*(log(N./Nc)+(N./Nc).*(64+0.05524.*(N./Nc).*(64+sqrt(N./Nc))).^(-1/4)); %quasi-Fermi level for CB
Efv=kB.*T.*(-1.*(log(P./Nv)+(P./Nv).*(64+0.05524.*(P./Nv).*(64+sqrt(P./Nv))).^(-1/4))-Eg); %quasi-Fermi level for VB (zero energy defined as CB minimum)
       
Eah=(Eg-E).*(me./(me+mhh)); %initial state in VB (heavy holes)
Eal=(Eg-E).*(me./(me+mlh)); %initial state in VB (light holes)
Ebh=(E-Eg).*(mhh./(me+mhh)); %final state in CB (heavy holes)
Ebl=(E-Eg).*(mlh./(me+mlh)); %final state in CB (light holes)
fc_Ebh=(1+exp((Ebh-Efc)./(kB.*T))).^(-1);%Fermi-Dirac probability of state Ebh being occupied by an electron
fc_Ebl=(1+exp((Ebl-Efc)./(kB.*T))).^(-1);%Fermi-Dirac probability of state Ebl being occupied by an electron
fv_Eah=(1+exp((Eah-Efv)./(kB.*T))).^(-1); %Fermi-Dirac probability of state Eah being occupied by an electron
fv_Eal=(1+exp((Eal-Efv)./(kB.*T))).^(-1); %Fermi-Dirac probability of state Eal being occupied by an electron
BF=(Chh./E).*sqrt(E-Eg).*(fv_Eah-fc_Ebh-1)+(Clh./E).*sqrt(E-Eg).*(fv_Eal-fc_Ebl-1); %change in alpha due to band filling effect (equation 11 in Bennet 1990)
% note, BF is in m^-1


end


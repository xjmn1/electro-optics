
function [Dn,Dk]=FK(Q,E,Eg,me,mhh,mlh,ri,C)
%Function to calculate change in real (Dn) and imaginary (Dk) parts of
%refractive index due to the Franz-Keldysh effect in semiconductors
%Based on the model provided in B. R. Bennett and R. A. Soref, “Electrorefraction and electroabsorption in InP, GaAs, GaSb, InAs, and InSb,” 
%IEEE Journal of Quantum Electronics, vol. 23, no. 12, pp. 2159–2166, Dec. 1987.

%Note, this function calls two other Matlab functions, F.m and G.m

%Input arguments:
%Q = Electric field (V m^-1)
%E = row vector of photon energies (J)
%Eg = material bandgap (J)
%me = electron effective mass (m0)
%mhh is heavy hole effective mass (m0)
%mlh = light hole effective mass (m0)
%ri is row vector of wavelength-dependent refractive index of material
%(same length as E) or scalar if fixed r.i. used
%C = fitted parameter obtained from fitting theory to experimental data
%(m^-1 J^1/2)(see Bennett, 1987)

% constants
c=2.998e8; % speed of light (m s^-1)
e=1.602e-19; % electron charge (C)
m0=9.109e-31; % electron mass (kg)
hbar=1.055e-34; % reduced Planck's constant (kg m^2 s^-1)

% convert variables to different units
omega=E./hbar; % convert E to angular frequency (s^-1)
omega_g=Eg./hbar; % angular frequency corresponding to band gap (s^-1)
mekg=me.*m0; % electron effective mass (kg)
mhhkg=mhh.*m0; % heavy hole effective mass (kg)
mlhkg=mlh.*m0; % light hole effective mass (kg)
mu_ehh=mekg.*mhhkg./(mekg+mhhkg); % electron-heavy hole reduced effective mass (kg)
mu_elh=mekg.*mlhkg./(mekg+mlhkg); % electron-light hole reduced effective mass (kg)
C=C./sqrt(hbar);


% calculate Dn and Dk

Theta_1=(e.^2.*Q.^2./(2.*mu_ehh.*hbar)).^(1/3); % equation 6a
Theta_2=(e.^2.*Q.^2./(2.*mu_elh.*hbar)).^(1/3); % equation 6b
       
x1=(omega_g-omega)./Theta_1; % equation 7
x2=(omega_g-omega)./Theta_2; % equation 7
B=ri.*c.*C./(mu_ehh.^(3/2).*(1+m0/mhh)+mu_elh.^(3/2).*(1+m0/mlh)); % B is calculated from eqn 28
Dn=real((B./(2.*ri.*omega.^2)).*(mu_ehh^(3/2).*(1+m0/mhh).*sqrt(Theta_1).*G(x1)+mu_elh.^(3/2).*(1+m0/mlh).*sqrt(Theta_2).*G(x2))); % equation 22
Dn(isnan(Dn))=0;
Dk=real((B./(2.*ri.*omega.^2)).*(mu_ehh^(3/2).*(1+m0/mhh).*sqrt(Theta_1).*F(x1)+mu_elh^(3/2).*(1+m0/mlh).*sqrt(Theta_2).*F(x2)));% equation 23
Dk(isnan(Dk))=0;        

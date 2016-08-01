% Script for calculating carrier-induced and Franz-Keldysh shifts in
% absorption and refractive index in semiconductors

%Based on the models for carrier induced effects provided in [1] and [2]
%and the model for the Franz-Keldysh effect provided in [3]

%% Variables:
%E = photon energy (J); row vector
%Q = row vector of electric field values (V m^-1)
%N = electron concentration (m^-3); row vector
%P = hole concentration (m^-3); row vector
%ri = real part of refractive index of material; scalar or row vector with
%same length as E
%me = electron effective mass (m0)
%mhh is heavy hole effective mass (m0)
%mlh = light hole effective mass (m0)
%Eg = material bandgap (J)
%C = fitted parameter from absorption data (m^-1 J^1/2)
%eps_s = static dielectic constant of material (units of vacuum permittivity eps_0)
%kappa = fitting parameter (J)
%type = 1 or 2: 1=n-doping or injected electron-hole plasma; 2=p-doping
%T = temperature (K)
%% Free carrier absorption

Eev=E./1.602e-19; %converts E from J into eV

DA_FCA=zeros(length(N),length(E));
DN_FCA=zeros(length(N),length(E));

for kk=1:length(N)

    [DA_FCA(kk,1:length(E)),DN_FCA(kk,1:length(E))]=FCabs(N(kk),P(kk),Eev,ri,me);
    % calculates change in absorption (m^-1) and r.i. due to free carrier
    % absorption
end

%% Bandgap shrinkage

BGS=BGshrink(me,eps_s,kappa,N,P,type);
%calculate bandgap shrinkage


DA_BGS=zeros(length(N),length(E));
DN_BGS=zeros(length(N),length(E));

for ii=1:length(N)
    
    DA_BGS(ii,1:length(E))=real((C./E).*sqrt(E-Eg-BGS(ii))-(C./E).*sqrt(E-Eg));
    %calculate change in absorption due to bandgap shrinkage using Eqn 19,
    %Bennett 1990
    DN_BGS(ii,1:length(E))=KKtrans(E,DA_BGS(ii,1:length(E)));
%calculate change in r.i. due to bandgap shrinkage using Kramers Kronig
%relation
    
end


%% Bandfilling

mekg=me.*9.109e-31; %converts me, mhh and mlh from m0 to kg
mhhkg=mhh.*9.109e-31; 
mlhkg=mlh.*9.109e-31;

Egshr=Eg+BGS; % convert Eg into row vector of bandgaps modified due to bandgap shrinkage
% row vector with same length as N

DA_BF=zeros(length(N),length(E));
DN_BF=zeros(length(N),length(E));

for jj=1:length(N)
    
    DA_BF(jj,1:length(E))=real(bandfilw(E,N(jj),P(jj),T,Egshr(jj),mekg,mhhkg,mlhkg,C));
    %calculate change in absorption due to bandfilling (using bandgaps modified to include shrinkage)
    
    DN_BF(jj,1:length(E))=KKtrans(E,DA_BF(jj,1:length(E)));
%calculate change in refractive index due to bandfilling
end



%% Combination of carrier effects

%add all three contributions together
DA_C=DA_FCA+DA_BGS+DA_BF;
DN_C=DN_FCA+DN_BGS+DN_BF;
DK_C=6.626e-34.*2.998e8.*DA_C./(4.*pi.*E); %convert absorption coefficient to imaginary part of r.i. 


%% Franz-Keldysh effect - can be run independently from the above sections

% Initialise arrays
DN_FK=zeros(length(Q),length(E));
DK_FK=zeros(length(Q),length(E));

for p=1:length(Q)
    [DN_FK(p,1:length(E)),DK_FK(p,1:length(E))]=FK(Q(1,p),E,Eg,me,mhh,mlh,ri,C);
end

%% Total shift
%note - only use this section if N,P and Q are all same length and each element of N,P and Q refer to the same layer/material 

DN=DN_C+DN_FK;
DK=DK_C+DK_FK;

%% References

%{
[1] B. R. Bennett, R. A. Soref, and J. A. Del Alamo, “Carrier-induced change in refractive index of InP, 
GaAs and InGaAsP,” IEEE Journal of Quantum Electronics, vol. 26, no. 1, pp.
113–122, Jan. 1990
[2] J.-P. Weber, “Optimization of the carrier-induced effective index change in InGaAsP waveguides-
application to tunable Bragg filters,” IEEE Journal of Quantum Electronics, vol. 30, no. 8, 
pp. 1801–1816, Aug. 1994.
[3] B. R. Bennett and R. A. Soref, “Electrorefraction and electroabsorption in InP, GaAs, GaSb, InAs, 
and InSb,” IEEE Journal of Quantum Electronics, vol. 23, no. 12, pp. 2159–2166, Dec. 1987.




%}

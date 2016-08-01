function Delta_n = KKtrans(E,Delta_alpha)
% Function to calculate change in real part of refractive index given change in optical absorption coefficient 
% using the Kramers-Kronig relation 
% Numerical integration is performed using Maclaurin's formula, see K. Ohta and H. Ishida, 
% “Comparison among Several Numerical Integration Methods for Kramers-Kronig Transformation,” 
% Appl Spectrosc, vol. 42, no. 6, pp. 952–957, Aug. 1988.

% Adapted from code provided by Mahmut Ruzi (c) 2015 via Mathworks File Exchange (see licence file text below)

% Input arguments:
% E and Delta_alpha can be row or column vectors but must have the same
% size and dimensions
% E is photon energy (J), Delta_alpha is change in absorption coefficient
% (m^-1)

c=2.998e8; 
hbar=1.055e-34; 


if size(E,2)>size(E,1)
    E=E'; 
end

if size(Delta_alpha,2)>size(Delta_alpha,1)
    Delta_alpha=Delta_alpha'; 
end

h=abs(E(1)-E(2)); % computation interval (points E(j) are assumed to be evenly spaced)

Delta_n=zeros(length(E),1);       % preallocate matrix for Delta_n (column vector)

for ii=1:2:length(E)               % odd indices
    Delta_n(ii)=(c.*hbar./pi).*2.*h.*sum(Delta_alpha(2:2:end)./(E(2:2:end).^2-(E(ii).*ones(length(E(2:2:end)),1)).^2));
end


for ii=2:2:length(E)                                % even indices
    Delta_n(ii)=(c.*hbar./pi).*2.*h.*sum(Delta_alpha(1:2:end)./(E(1:2:end).^2-(E(ii).*ones(length(E(1:2:end)),1)).^2));
end


end

%Licence file  text:
%{
Copyright (c) 2015, Mahmut Ruzi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
%}

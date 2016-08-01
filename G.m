function [G] = G(x)
% Electrooptic function G(x) as defined in equation 5, B. R. Bennett and R. A. Soref, 
% “Electrorefraction and electroabsorption in InP, GaAs, GaSb, InAs, and InSb,” 
% IEEE Journal of Quantum Electronics, vol. 23, no. 12, pp. 2159–2166, Dec. 1987.

if x<0
    H=0;
else
    H=1;
end

G=pi().*(airy(1,x).*airy(3,x)-x.*airy(x).*airy(2,x))+sqrt(x).*H;

end


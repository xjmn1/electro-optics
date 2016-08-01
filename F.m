function [F] = F(x)
% Electrooptic function F(x) as defined in equation 5, B. R. Bennett and R. A. Soref, 
% “Electrorefraction and electroabsorption in InP, GaAs, GaSb, InAs, and InSb,” 
% IEEE Journal of Quantum Electronics, vol. 23, no. 12, pp. 2159–2166, Dec. 1987.

if -x<0
    H=0;
else
    H=1;
end

F=pi().*((abs(airy(1,x))).^2-x.*(airy(x)).^2)-sqrt(-x).*H;

end


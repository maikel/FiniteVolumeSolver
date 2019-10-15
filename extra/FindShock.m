function [max_drho, maxx, maxi] = FindShock(X, rho)
   drho = abs(rho(2:end) - rho(1:end-1));
   [max_drho, maxi] = max(drho);
   maxx = X(maxi);
end


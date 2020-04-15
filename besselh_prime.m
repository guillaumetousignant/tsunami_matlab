function [result] = besselh_prime(nu, z)
    result = nu * besselh(nu, z)/z - besselh(nu+1, z);
end
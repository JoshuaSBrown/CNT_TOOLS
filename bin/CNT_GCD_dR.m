function [ dR ] = CNT_GCD_dR( n, m )
%Calculates greatest common devisor of (2n+m) and (2m+n) given the chiral
%vectors of a CNT ( n, m)
%   Detailed explanation goes here

%Checks
assert(0<=n,'n must be greater than 0');
assert(m<=n,'m must be <= n');
assert(m>=0,'m must be >= 0');
assert(rem(n,1)==0,'n must be a whole number');
assert(rem(m,1)==0,'m must be a whole number');

%Calculation
dR = gcd((2*n+m), (2*m+n) );

end


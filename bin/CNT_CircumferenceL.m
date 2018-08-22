function [ L ] = CNT_CircumferenceL( n, m )
%Given a CNT's chiral vectors n & m this function will calculate the
%circumferential Length L
%   To call the function:
%   [ L ] = CNT_CircumferenceL( 4, 3 )
%   L is the circumferenctial Length and is returned as a scaler
%   n = 4 and m = 3 in the case above and represents the Chiral vectors

%Checks
assert(0<=n,'n must be greater than 0');
assert(m<=n,'m must be <= n');
assert(m>=0,'m must be >= 0');
assert(rem(n,1)==0,'n must be a whole number');
assert(rem(m,1)==0,'m must be a whole number');

%Constants
a=1.42*(3)^(1/2);

%Calculation
L = a*(n^2+m^2+n*m)^(1/2);

end


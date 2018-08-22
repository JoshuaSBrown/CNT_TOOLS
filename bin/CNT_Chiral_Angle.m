function [ theta ] = CNT_Chiral_Angle( n, m )
%Calculates the angle between the Chiral vector Cvec and the unit vector a1
%where: 
% Cvec = na1 + ma2
% a1 = [ (3)^(1/2)/2 a,  a/2 ]
% a2 = [ (3)^(1/2)/2 a,  -a/2 ]
%   theta should be returned in degrees and range between 0 and 30: 
%   0 - zigzag
%   30 - armchair
%   0 < theta < 30 - chiral

%Checks
assert(0<=n,'n must be greater than 0');
assert(m<=n,'m must be <= n');
assert(m>=0,'m must be >= 0');
assert(rem(n,1)==0,'n must be a whole number');
assert(rem(m,1)==0,'m must be a whole number');

%Calculation
theta = acosd((2*n+m)/(2*(n^2+m^2+n*m)^(1/2)));

%Post Check
assert(theta>=0,'theta must be greater or equal to 0');
assert(theta<=30.01,'theta must be less than or equal to 30');

end


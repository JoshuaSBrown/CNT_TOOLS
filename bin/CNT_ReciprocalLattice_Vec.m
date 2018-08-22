function [ K1, K2 ] = CNT_ReciprocalLattice_Vec( n, m )
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here


%Checks
assert(0<=n,'n must be greater than 0');
assert(m<=n,'m must be <= n');
assert(m>=0,'m must be >= 0');
assert(rem(n,1)==0,'n must be a whole number');
assert(rem(m,1)==0,'m must be a whole number');

%Constants
a=1.42*(3)^(1/2);

%Vectors
b1=[ 2*pi/a,  2*pi/((3)^(1/2)*a)];
b2=[ -2*pi/a, 2*pi/((3)^(1/2)*a)];

%Calculations
[ N ] = CNT_UnitCell_Num_Hex( n, m);

[ t1, t2] = CNT_Translational_Vec_t1t2(n,m);

K1 = 1/N*[ -t2*b1(1)+t1*b2(1), -t2*b1(2)+t1*b2(2)];
K2 = 1/N*[ m*b1(1)-n*b2(1), m*b1(2)-n*b2(2)];

end


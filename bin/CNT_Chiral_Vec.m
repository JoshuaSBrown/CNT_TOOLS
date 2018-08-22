function [ Cvec ] = CNT_Chiral_Vec( n, m )
%This function takes the chiral vectors n & m and returns the chiral vector
%Cvec
%   To run this function simply type:
%   [ Cvec ]= ChiralVec( 4, 3 ); 
%   Where n = 4 and m = 3;

%Checks
assert(0<n,'n must be greater than 0');
assert(m<=n,'m must be <= n');
assert(m>=0,'m must be >= 0');
assert(rem(n,1)==0,'n must be a whole number');
assert(rem(m,1)==0,'m must be a whole number');

%Constants
a=1.42*(3)^(1/2);

%Vectors
a1=[a/2, (3)^(1/2)/2*a];
a2=[-a/2, (3)^(1/2)/2*a];

%Calculation 
%Cvec = n*a1 + m*a2
Cvec=[n*a1(1)+m*a2(1), n*a1(2)+m*a2(2)];

end


function [ Rvec ] = CNT_Symmetry_Vec( n, m )
%Calculates the Symmetry vector Rrvec. It is the site vector having the
%smallest component in the direction of the Chiral vector
%   Detailed explanation goes here

%Constants
a=1.42*(3)^(1/2);

%Vectors
a1=[(3)^(1/2)/2*a, a/2];
a2=[(3)^(1/2)/2*a, -a/2];

%Calculations
[ p, q] = CNT_Symmetry_Vec_pq( n, m);

Rvec = [p*a1(1)+p*a2(1),q*a1(2)+q*a2(2)];

end


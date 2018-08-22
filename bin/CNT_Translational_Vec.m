function [ T ] = CNT_Translational_Vec( n, m )
%Calculates the translational vector of a CNT given (n, m) chiral vectors
%the Translational Vector is a unit vector of the CNT and runs parallel to 
%the axis of the CNT and perpendicular to the Chiral Vector 
%   Detailed explanation goes here

%Constants
a=1.42*(3)^(1/2);

%Vectors
a1=[(3)^(1/2)/2*a, a/2];
a2=[(3)^(1/2)/2*a, -a/2];

%Calculations
[ t1, t2 ] = CNT_Translational_Vec_t1t2( n, m);

T = [ t1*a1(1)+t2*a2(1), t1*a1(2)+t2*a2(2) ];

end


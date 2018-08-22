function [ t1, t2 ] = CNT_Translational_Vec_t1t2( n, m )
%Calculates t1 and t2 given the chiral vectors ( n, m) of a CNT. The
%Translational vector may be expressed in terms of the basis vectors a1 and
%a2 by multiplying them by the appropriate integers t1 and t2 respectively
%   Detailed explanation goes here

%Calculations
[ dR ] = CNT_GCD_dR(n,m);

t1 = (2*m+n)/dR;
t2 = -(2*n+m)/dR;

%Post Check
assert(rem(t1,1)==0,'t1 has been calculated and is not a whole number');
assert(rem(t2,1)==0,'t2 has been calculated and is not a whole number');

end


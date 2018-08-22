function [ N ] = CNT_UnitCell_Num_Hex( n, m)
%Calculates the number of hexagons (defined by the lattice structure of
%graphene) will fit in a unit cell of a CNT which is defined by the
%Translational Vector Tvec and Chiral Vector Cvec.
%   Detailed explanation goes here

%Constants
a = 1.42*(3)^(1/2);

%Calculations
[ L ] = CNT_CircumferenceL( n, m);

[ dR ] = CNT_GCD_dR( n, m);

N = 2*L^2/(a^2*dR);

end


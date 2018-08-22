function [ Num_C_Atoms ] = CNT_UnitCell_Num_Atoms( n, m )
%Calculates the number of Carbon atoms in a unit cell of a CNT. The CNT
%unit cell is defined by the Translational Vector Tvec and the Chiral
%Vector Cvec. There are two carbon atoms for every hexagon in a CNT unit
%cell. 
%   Detailed explanation goes here

%Calculations
[ N ] = CNT_UnitCell_Num_Hex( n, m);

Num_C_Atoms = 2*N;

end


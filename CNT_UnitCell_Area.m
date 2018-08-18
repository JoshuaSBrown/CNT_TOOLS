function [ Area ] = CNT_UnitCell_Area( n, m )
%Calculates the unit cell of a (n,m) CNT
%   Area should be returned in units of Ang^2

[ Cvec ] = CNT_Chiral_Vec( n, m );
vec1 = [ Cvec(1) Cvec(2) 0];

[ T ] = CNT_Translational_Vec( n, m );
vec2 = [ T(1) T(2) 0];

vec3 = cross(vec1, vec2);

Area = ( sum(vec3.^2) )^(1/2);
end


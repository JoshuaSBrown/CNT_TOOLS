function [ theta ] = CNT_Rotation_Angle( n, m )
%The rotation angle describes the relationship between the two
%Coordinate systems (kx,ky) and (kT, kC)

theta = abs(atan(-(n-m)/((3)^(1/2)*(n+m))));

end


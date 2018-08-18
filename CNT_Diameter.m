function [ CNT_D ] = CNT_Diameter( n, m )
%Given the chiral vectors this function calculates the CNT diameter. 
%   This function works by calling the CNT_Circumference function and then
%   dividing by pi

%Checks
[ L ] = CNT_CircumferenceL( n, m);

%Calculation
CNT_D = L/pi;

end


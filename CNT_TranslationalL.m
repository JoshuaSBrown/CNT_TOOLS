function [ TL ] = CNT_TranslationalL( n, m )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[T]=CNT_Translational_Vec(n,m);

TL = (T(1)^2+T(2)^2)^(1/2);

end


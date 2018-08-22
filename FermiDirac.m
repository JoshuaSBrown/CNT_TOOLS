function [ F ] = FermiDirac( DeltaE, T )
%Fermi Dirac distribution
%   DeltaE is the E-Ef where Ef is the fermi energy of
%   the material in [eV] DeltaE can be submitted as a 
%   matrix or a single element and T is the temperature

Kb = 8.617332478*10^-5; % [eV/K]

F = 1./(exp(DeltaE/(Kb*T))+1);

end


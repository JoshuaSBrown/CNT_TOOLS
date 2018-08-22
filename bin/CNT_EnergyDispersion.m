function [ h ] = CNT_EnergyDispersion( n, m )
%Plot's the Energy dispersion curves for a nanotube
%   Detailed explanation goes here

%Checks
assert(0<n,'n must be greater than 0');
assert(m<=n,'m must be <= n');
assert(m>=0,'m must be >= 0');
assert(rem(n,1)==0,'n must be a whole number');
assert(rem(m,1)==0,'m must be a whole number');

%Constants
a=1.42*(3)^(1/2);
Resolution=1000;
%t is the nearest neighbors overlap energy nominally between 2.5 and 3.2
%[eV] This value has been reported by J. W. G. Wilder and C.
%Dekker et al. Nature (1998) to be approximately 2.7+-0.1 [eV]
t= 2.7;
[ T ] = CNT_Translational_Vec( 5,5);
Boundary=pi/(T(1)^2+T(2)^2)^(1/2);
[ N ] = CNT_UnitCell_Num_Hex( n, m);
[ t1, t2 ] = CNT_Translational_Vec_t1t2( n, m);

%Vectors
b1 = 2*pi/a*[1/(3)^(1/2), 1 ];
b2 = 2*pi/a*[1/(3)^(1/2), -1 ];

b1(1)=(3)^2/2*b1(1)-1/2*b1(2);
b1(2)=1/2*b1(1)+(3)^(1/2)/2*b1(2);
b2(1)=(3)^2/2*b2(1)-1/2*b2(2);
b2(2)=1/2*b2(1)+(3)^(1/2)/2*b2(2);
%Initialization
E=zeros(Resolution,2*int64(N));
k=linspace(-Boundary,0,Resolution);

%Calculation

K1(1)=1/N*(-t2*b1(1)+t1*b2(1));
K1(2)=1/N*(-t2*b1(2)+t1*b2(2));

K2(1)=1/N*(m*b1(1)-n*b2(1));
K2(2)=1/N*(m*b1(2)-n*b2(2));

K2Ab=(K2(1)^2+K2(2)^2)^(1/2);
%Armchair
% if (n==m)
%    for q=1:2*n
%       E(:,q)=t*(1+4*cos(q*pi/n)*cos(k*a/2)+4*cos(k*a/2).^2).^(1/2); 
%       E(:,q+2*n)=-t*(1-4*cos(q*pi/n)*cos(k*a/2)+4*cos(k*a/2).^2).^(1/2);
%    end
% %Zigzag
% elseif(m==0)
%    for q=1:2*n
%       E(:,q)=t*(1+4*cos(q*pi/n)*cos((3)^(1/2)*k*a/2)+4*cos(q*pi/n).^2).^(1/2); 
%       E(:,q+2*n)=-t*(1-4*cos(q*pi/n)*cos((3)^(1/2)*k*a/2)+4*cos(q*pi/n).^2).^(1/2);
%    end
% %Chiral
% else
    for q=1:N
        
       %kx=(m-n)*k./(2*(m^2+n^2+m*n)^(1/2))+q*2*pi/(a*N)*(t1-t2)/(3)^(1/2);
       %ky=(3)^(1/2)*k./2*(m+n)/(m^2+n^2+m*n)^(1/2)-q*2*pi/(a*N)*(t1+t2)/(3)^(1/2);
       kx=k.*K2(1)/K2Ab+q*K1(1);
       ky=k.*K2(2)/K2Ab+q*K1(2);
       E(:,q)=t*(1+4*cos((3)^(1/2)*kx*a/2).*cos(ky*a/2)+4*cos(ky*a/2).^2).^(1/2);
       E(:,q+int64(N))=-t*(1+4*cos((3)^(1/2)*kx*a/2).*cos(ky*a/2)+4*cos(ky*a/2).^2).^(1/2);
    end
% end

%Plotting
h=figure(1);
plot(k,E);
xlim([-Boundary 0]);

end


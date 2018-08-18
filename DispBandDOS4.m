function [carrierConc, carrierConc2, CNT_D ] = DispBandDOS4( n, m)
%CNT_D returned in units of Ang
%carrierConc and carrierConc2 are both in units of 1/cm3
%carrierConc - for a layer of the CNT, (Might need the layer if considering
%MWCNT)
%carrierConc2 - for a whole SWCNT
close all

Temp = 300;      % [K]
Elimit=-1;       %Only consider dispersion relation where bands cross below Elimit
%If set equal to -1 all bands are included
Resolution=4000;

a=1.42*(3)^(1/2); %Units of Angstroms

a1=[a/2,(3)^(1/2)*a/2];         % a is in units of [Ang]
a2=[ -a/2,(3)^(1/2)*a/2];

Ch=[n*a1(1)+m*a2(1),n*a1(2)+m*a2(2)];   % [Ang] This has been verified with Dressler

[ dR ] = CNT_GCD_dR(n,m);

t1 = (2*m+n)/dR;
t2 = -(2*n+m)/dR;

T = [ t1*a1(1)+t2*a2(1), t1*a1(2)+t2*a2(2) ]; %[Ang]

L=a*sqrt(n^2+m^2+n*m);                  %Angstroms
Tmax=(3)^(1/2)*L/dR;                    %Angstroms
%From Dresselhaus pg 40 verified

[ N ] = CNT_UnitCell_Num_Hex( n, m);
%From Dresselhaus verified with table 3.2
q=1:round(N);

figure(3);
x=linspace(0,a1(1),20);
y=linspace(0,a1(2),20);
hold on
plot(x,y,'r');
x=linspace(0,a2(1),20);
y=linspace(0,a2(2),20);
plot(x,y,'r');
x=linspace(0,Ch(1),100);
y=linspace(0,Ch(2),100);
plot(x,y);
x=linspace(0,T(1),100);
y=linspace(0,T(2),100);
plot(x,y);
line([a/2 a/2],[-1.42/2, 1.42/2],'Color','bl');
line([-a/2 -a/2],[-1.42/2, 1.42/2],'Color','bl');
line([-a/2 0],[1.42/2 (3)^(1/2)*a/2-1.42/2],'Color','bl')
line([ 0 a/2],[(3)^(1/2)*a/2-1.42/2 1.42/2],'Color','bl')
axis equal
set(gcf,'Color','white');
%Jorgenson Fig 4 pg 9

b1=[ 2*pi/(sqrt(3)*a),  2*pi/(a)];
b2=[ 2*pi/((3)^(1/2)*a), -2*pi/a];
K1=((2*n+m)*b1+(2*m+n)*b2)/(N*dR);   %[1/Ang]
K2=(m*b1-n*b2)/N;   %[1/Ang]
%From Marulanda

b1=[ 2*pi/a,  2*pi/((3)^(1/2)*a)];
b2=[ -2*pi/a, 2*pi/((3)^(1/2)*a)];
%Jorgenson pg 9

V2=Ch(1)*T(2)-T(1)*Ch(2);           % [Ang^2]
GT=[2*pi*T(2)/V2, -2*pi*T(1)/V2];   % [1/Ang]
Gc=[-2*pi*Ch(2)/V2, 2*pi*Ch(1)/V2];

figure(4);
x=linspace(0,b1(1),20);
y=linspace(0,b1(2),20);
hold on
plot(x,y,'r');
x=linspace(0,b2(1),20);
y=linspace(0,b2(2),20);
plot(x,y,'r');
x=linspace(0,Gc(1),100);
y=linspace(0,Gc(2),100);
plot(x,y);
x=linspace(0,GT(1),100);
y=linspace(0,GT(2),100);
plot(x,y);
axis equal
set(gcf,'Color','white');
%Jorgenson Fig 4 pg 9
knew=linspace(-pi/Tmax,pi/Tmax,(Resolution*2+1));     %[1/Ang]
Knew=zeros(Resolution*2+1,round(N),2);
for i=1:round(N)
    Knew(:,i,1)=knew.*K2(1)/sqrt(K2(1)^2+K2(2)^2)+q(i)*K1(1);  %[1/Ang]
    Knew(:,i,2)=knew.*K2(2)/sqrt(K2(1)^2+K2(2)^2)+q(i)*K1(2);
    %From Marulanda
end

%From Marulanda
figure(9);
hold on;
xlabel('k');
ylabel('Energy');
Eini=zeros((Resolution*2+1),round(N));
Vpp=2.5;        %Units of eV
%Teri W. Odom

dEdk=zeros(1,(Resolution+1));
DOS=zeros(Resolution,1);

Nnew=0;
for i=1:round(N)
    Eini(:,i)=Vpp*(1+4*cos(sqrt(3)*Knew(:,i,1)*a/2).*cos(Knew(:,i,2)*a/2)+4*cos(Knew(:,i,2)*a/2).^2).^(1/2);  %Units [eV]
    if(min(Eini(:,i))<Elimit)
        Nnew=Nnew+1;
    end
    plot(knew,Eini(:,i));
end

set(gcf,'color','w');

%Check to see if Dispersion Folder Exists
if (exist('Dispersion_Library','dir')~=7)
    mkdir('Dispersion_Library');
end 

str2 = strcat('Dispersion_Library/Disp_',num2str(n),'_',num2str(m),'.fig');
if (exist(str2,'file')~=2)
    saveas(gcf,str2);
end

figure(10);
hold on;
xlabel('k');
ylabel('Energy');

if (Elimit==-1)
    Enew=Eini;
else
    j=1;
    Enew=zeros((Resolution*2+1),Nnew);
    for i=1:N
        if(min(Eini(:,i))<Elimit)
            Enew(:,j)=Eini(:,i);
            plot(knew,Enew(:,j));
            j=j+1;
            
        end
    end
    
    N=Nnew;
end

Emax=max(max(Enew));

figure(11);
hold on;
xlabel('Energy [eV]');
ylabel('DOS');

%The DOS is the length of Resolution
%Starting at Energy 0 and going to Emax
inc=(Emax)/Resolution;

i=1;
p=1;
while i<=round(N)
    
    EnewUpper = Enew(Resolution:Resolution*2+1,i);
    
    
    if(min(EnewUpper(1:Resolution)-EnewUpper(2:Resolution+1))<0 && max(EnewUpper(1:Resolution)-EnewUpper(2:Resolution+1))>0 )
        %Finding the element where the switch occurs
        value = EnewUpper(1)-EnewUpper(2);
        if (value>0)
            Default=1;
        else
            Default=-1;
        end
        
        Switch = Default;
        
        j=2;
        m=1;
        store=zeros(1,1);
        
        while (j<Resolution)
            
            value = EnewUpper(j)-EnewUpper(j+1);
            
            if (value>0)
                Switch=1;
            elseif(value<0)
                Switch=-1;
            end
            %If value is 0 Switch does not change
            
            %Store switching points
            if (Switch~=Default && m==1)
                if(j==2)
                    Default=-Default;
                else
                    store(1)=j-1;
                    m=m+1;
                    Default = Switch;
                end
            elseif (Switch~=Default && j~=2)
                store = [store; j-1];
                m = m+1;
                Default = Switch;
            end
            j=j+1;
        end
        
        m=1;
        Ini=1;
        
        if (length(store)>=1 && store(1)~=0)
            while(m<=length(store))
                dEdk2 = (EnewUpper(Ini+1:store(m))-EnewUpper(Ini:store(m)-1))'./(knew(Ini+1:store(m))-knew(Ini:store(m)-1));       %[eV Ang]
                Emid=(EnewUpper(Ini+1:store(m))+EnewUpper(Ini:store(m)-1))/2;
                Emin1=min(Emid);
                Emax1=max(Emid);
                
                %If there are not enough data points in the cut section
                %just ignore it. This will only introduce a negligable
                %amoutn of error. 
                Ni=ceil(Emin1/inc);
                Nf=floor(Emax1/inc);
                if((Nf-Ni)/inc>7)
                    [vq, Ni, Nf ] = FitDOS(Emin1, Emax1, inc, Emid, dEdk2);
                    DOS(Ni:Nf,1)=DOS(Ni:Nf,1)+vq';  % [1/( eV Ang )]
                end
                Ini = store(m)+1;
                m=m+1;
                p=p+1;
            end
            
            dEdk2 = (EnewUpper(Ini+1:Resolution+1)-EnewUpper(Ini:Resolution))'./(knew(Ini+1:Resolution+1)-knew(Ini:Resolution));
            Emid=(EnewUpper(Ini+1:Resolution+1)+EnewUpper(Ini:Resolution))/2;
            Emin1=min(Emid);
            Emax1=max(Emid);
            Ni=ceil(Emin1/inc);
            Nf=floor(Emax1/inc);
            if((Nf-Ni)/inc>7)
                [vq, Ni, Nf ] = FitDOS(Emin1, Emax1, inc, Emid, dEdk2);
                DOS(Ni:Nf,1)=DOS(Ni:Nf,1)+vq'; % [1/( eV Ang )]
            end
        else
            dEdk(1,:)=(Enew((Resolution+1):(Resolution*2+1),i)-Enew(Resolution:2*Resolution,i))'./(knew(Resolution+1:Resolution*2+1)-knew(Resolution:2*Resolution));
            Emid=(Enew((Resolution+1):(Resolution*2+1),i)+Enew(Resolution:2*Resolution,i))/2;
            Emin1=min(Emid);
            Emax1=max(Emid);
            [vq, Ni, Nf ] = FitDOS(Emin1, Emax1, inc, Emid, dEdk);
            DOS(Ni:Nf,1)=DOS(Ni:Nf,1)+vq';   % [1/( eV Ang )]
        end
    else
        
        dEdk(1,:)=(Enew((Resolution+1):(Resolution*2+1),i)-Enew(Resolution:2*Resolution,i))'./(knew(Resolution+1:Resolution*2+1)-knew(Resolution:2*Resolution));
        Emid=(Enew((Resolution+1):(Resolution*2+1),i)+Enew(Resolution:2*Resolution,i))/2;
        Emin1=min(Emid);
        Emax1=max(Emid);
        [vq, Ni, Nf ] = FitDOS(Emin1, Emax1, inc, Emid, dEdk);
        DOS(Ni:Nf,1)=DOS(Ni:Nf,1)+vq'; % [1/( eV Ang )]
    end
    
    i=i+1;
    p=p+1;
end
DOS=2*2*DOS/(N);    % DOS [ 1 / (eV Ang) ]  Take out GTmax [Ang] = [1/eV Ang]
%The first 2 is for the spin of the electrons thus the DOS is doubled for
%each band, the second 2 is for the -k values we only integrated over the
%positive side. 
%Plot Graphene
%x=linspace(0,2*a,Resolution);
figure;
hold on
plot(inc/2:inc:Emax-inc/2,DOS)
xlabel('Energy [eV]');
ylabel('[1/eV]');
if (Elimit==-1)
    xlim([0 Emax]);
else
    xlim([0 Elimit]);
end

%print DOS to file
printDOS( n, m, inc/2:inc:Emax-inc/2, DOS);

[ F ] = FermiDirac( inc/2:inc:Emax-inc/2, Temp );

figure;
hold on
plot(inc/2:inc:Emax-inc/2,F);
xlabel('[eV]');
ylabel('Occupied States');
figure;
hold on
plot(inc/2:inc:Emax-inc/2,DOS.*F');
[ CNT_D ] = CNT_Diameter( n, m );                                 %Units [Ang]
Area = (CNT_D/2*10^-8)^2*pi-((CNT_D/2-0.617)*10^-8)^2*pi;         %Units [cm^2]
Area2 = (CNT_D/2*10^-8)^2*pi;                                     %Units [cm^2]
%This is the carrier concentration using the Area that Marulanda uses
carrierConc = trapz(inc/2:inc:Emax-inc/2,DOS.*F'/(Area*10^-8));   %[1/cm^3]

carrierConc2 = trapz(inc/2:inc:Emax-inc/2,DOS.*F'/(Area2*10^-8));

display('Carrier concentration in units of [1/cm^3]');
display(carrierConc)
end

function [vq, Ni, Nf ] = FitDOS(Emin1, Emax1, inc, Emid, dEdk)

Ni=ceil(Emin1/inc);
Nf=floor(Emax1/inc);
intervals=Ni*inc:inc:Nf*inc;

lx=find(Emid==min(Emid));
[rx,~]=ind2sub(size(Emid),lx);
Bot = rx(1);

lx=find(Emid==max(Emid));
[rx,~]=ind2sub(size(Emid),lx);
Top = rx(1);

Temp = sort(Emid);
MiddleNum = Temp(round(length(Emid)/2));
lx=find(Emid==MiddleNum);
[rx,~]=ind2sub(size(Emid),lx);
Mid = rx(1);
slope=(abs(1/dEdk(1,Top))-abs(1/dEdk(1,Bot)))/(Emid(Top)-Emid(Bot));
slope2=(abs(1/dEdk(1,Top))-abs(1/dEdk(1,Mid)))/(Emid(Top)-Emid(Mid));
slope3=(abs(1/dEdk(1,Mid))-abs(1/dEdk(1,Bot)))/(Emid(Mid)-Emid(Bot));

rval1 = 0;
rval2 = 0;
rval3 = 0;
rval4 = 0;
if(slope2>0 && slope3<0 )
    
    [ vq, h, rval1, rval2, rval3 ] = fitMiddle( Emid',abs(1./dEdk), intervals);
    legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
    set(gcf,'color','w')
    figure(11)
    hold on
    plot(intervals,vq);
end

if(rval1<0.93 && rval2<0.93 && rval3<0.93 && rval4<0.93)
    if (slope>0)
       
        [ vq, h ] = fitRight(Emid',abs(1./dEdk), intervals);
        legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
        figure(11)
        hold on
        plot(intervals,vq);
    else
        
        [ vq, h ] = fitLeft( x, y, intervals);
        legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
        figure(11)
        hold on
        plot(intervals,vq);
       
    end
end

end

function [ vq, h, rval1, rval2, rval3, rval4 ] = fitMiddle( x, y, intervals)

[xData, yData] = prepareCurveData( x, y );

%Fit 1
% Set up fittype and options.
ft = fittype( 'a/(abs(b-x)^n*abs((c-x))^n)', 'independent', 'x', 'dependent', 'y' );
opts1 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts1.Display = 'Off';
opts1.Lower = [-Inf min(x)*0.9 max(x)*1.00000000001 -Inf];
opts1.Robust = 'Bisquare';
opts1.StartPoint = [0.860114455139906 min(x)*0.99 max(x)*1.01 0.141886338627215];
opts1.Upper = [Inf min(x)*0.999999999 max(x)*1.1 Inf];
% Fit model to data.
[fitresult1, gof] = fit( xData, yData, ft, opts1 );
rval1 = gof.rsquare;
Coeffs1 = coeffvalues(fitresult1);
  
%Fit 2
% Set up fittype and options.
ft = fittype( 'a*x/(abs(x^2-c^2)*abs(x^2-b^2))^(1/2)', 'independent', 'x', 'dependent', 'y' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.MaxIter = 800;
opts2.Robust = 'Bisquare';
opts2.Lower = [-Inf min(x)*0.99 max(x)*1.000000001];
opts2.Upper = [Inf min(x)*0.999999 max(x)*1.03];
rval2=0;
inc=0;
while (rval2<0.93 && inc<5)
    opts2.StartPoint = [rand min(x)*0.999 max(x)*1.01];
    % Fit model to data.
    [fitresult2, gof] = fit( xData, yData, ft, opts2 );
    inc=inc+1;
    rval2 = gof.rsquare;
end
Coeffs2 = coeffvalues(fitresult2);

%Fit 3
 % Set up fittype and options.
ft = fittype( 'a*x/(abs(x^n-b^n)*abs(x^m-c^m))^(1/2)', 'independent', 'x', 'dependent', 'y' );
opts3 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts3.Display = 'Off';
opts3.Lower = [-Inf min(x)*0.9 max(x)*1.000001 1 1];
opts3.Robust = 'Bisquare';
opts3.Upper = [Inf min(x)*0.9999999 max(x)*1.1 Inf Inf];
inc =1;
rval3 = 0;
while (inc<5 && rval3<0.93)
    opts3.StartPoint = [rand min(x)*0.99 max(x)*1.01 2 2];
    % Fit model to data.
    [fitresult3, gof] = fit( xData, yData, ft, opts3 );
    rval3 = gof.rsquare;
    inc=inc+1;
end
Coeffs3 = coeffvalues(fitresult3);

%Fit 4
% Set up fittype and options.
ft = fittype( 'poly5' );
% Fit model to data.
[fitresult4, gof] = fit( xData, yData, ft );
Coeffs4 = coeffvalues(fitresult4);
rval4 = gof.rsquare;

if(rval1>=0.93 && rval1>=rval2 && rval1>=rval3 && rval1>=rval4)
    vq = abs(Coeffs1(1)./(abs(Coeffs1(2)-intervals).^Coeffs1(4).*abs(Coeffs1(3)-intervals).^Coeffs1(4)));
    %figure;
    %h = plot( fitresult1, x, y );
elseif(rval2>=0.93 && rval2>=rval3 && rval2>=rval4)
    vq = Coeffs2(1)*intervals./(abs(intervals.^2-Coeffs2(2)^2).*abs(intervals.^2-Coeffs2(3)^2)).^(1/2);%             vq = Coeffs3(2)+exp(-Coeffs3(1)*(intervals-Coeffs3(3)));
    %figure;
    %h = plot( fitresult2, x,y );
elseif(rval3>0.93 && rval3>=rval4)
    vq = Coeffs3(1)*intervals./(abs(intervals.^Coeffs3(5)-Coeffs3(2)^Coeffs3(5)).*abs(intervals.^Coeffs3(4)-Coeffs3(3)^Coeffs3(4))).^(1/2);%             vq = Coeffs3(2)+exp(-Coeffs3(1)*(intervals-Coeffs3(3)));
    %figure;
    %h = plot( fitresult3, x, y );
elseif( rval4>0.93)
    vq = Coeffs4(1)*(intervals).^5+Coeffs4(2)*(intervals).^4+Coeffs4(3)*(intervals).^3+Coeffs4(4)*(intervals).^2+Coeffs4(5)*intervals+Coeffs4(6);
    %figure;
    %h = plot( fitresult4, x, y );
end
end

function [vq, h ] = fitRight(x, y, intervals)

[xData, yData] = prepareCurveData( x, y );

%Fit attempt 1
rval2=0;
if(length(x)>6)
    % Set up fittype and options.
    ft = fittype( 'poly5' );
    % Fit model to data.
    [fitresult2, gof] = fit( xData, yData, ft );
    Coeffs2 = coeffvalues(fitresult2);
    rval2 = gof.rsquare;
end

%Fit attempt 2
str=strcat('a*x/(abs(',num2str(max(x)*1.00000001,12),'^2-x^2))^(1/2)');
% Set up fittype and options.
ft = fittype( str, 'independent', 'x', 'dependent', 'y' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Algorithm = 'Levenberg-Marquardt';
opts2.Display = 'Off';
opts2.Robust = 'Bisquare';
inc=0;
rval3=0;
while (rval3<0.9 && inc<5)
    opts2.StartPoint = rand;
    % Fit model to data.
    [fitresult3, gof] = fit( xData, yData, ft, opts2 );
    inc=inc+1;
    rval3=gof.rsquare;
end
Coeffs3 = coeffvalues(fitresult3);

%Fit 3
str=strcat('a*x/(abs(',num2str(max(x)*1.00000001,12),'^2-x^2))^(1/2)');
% Set up fittype and options.
ft = fittype( str, 'independent', 'x', 'dependent', 'y' );
opts3 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts3.Display = 'Off';
opts3.Lower = 0.0001 ;
opts3.Robust = 'Bisquare';
opts3.Upper = 1000000 ;
rval4=0;
inc=0;
while(rval4<0.9 && inc<5)
    opts3.StartPoint = 1000*rand+0.0001;
    % Fit model to data.
    [fitresult4, gof] = fit( xData, yData, ft, opts3 );
    rval4 = gof.rsquare;
    inc=inc+1;
end
Coeffs4 = coeffvalues(fitresult4);

%Fit 4
% Set up fittype and options.
ft = fittype( 'a*x/(abs(b^2-x^2))^(1/2)', 'independent', 'x', 'dependent', 'y' );
opts4 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts4.Display = 'Off';
bstart = max(x);
opts4.Lower = [0.00001 bstart*1.00001];
opts4.Robust = 'Bisquare';
opts4.Upper = [10000 bstart*1.5];
rval5=0;
inc=0;
while(rval5<0.9 && inc<5)
    opts4.StartPoint = [10*rand+0.01 bstart*1.3];
    % Fit model to data.
    [fitresult5, gof] = fit( xData, yData, ft, opts4 );
    rval5 = gof.rsquare;
    inc=inc+1;
end
Coeffs5 = coeffvalues(fitresult5);

if ( rval2>=rval3 && rval2>=rval4 && rval2>=rval5 )
    vq = Coeffs2(6)+Coeffs2(5)*intervals+Coeffs2(4)*intervals.^2+Coeffs2(3)*intervals.^3+Coeffs2(2)*intervals.^4+Coeffs2(1)*intervals.^5;
    %figure;
    %h = plot( fitresult2, x, y );
elseif ( rval3>=rval2 && rval3>=rval4 && rval3>=rval5 )
    vq = Coeffs3(1)*intervals./(abs((max(x)*1.00000001)^2-intervals.^2)).^(1/2);
    %figure;
    %h = plot( fitresult3, x, y );
elseif ( rval4>=rval2 && rval4>=rval3 && rval4>=rval5 )
    vq = Coeffs4(1)*intervals./(abs((max(x)*1.00000001)^2-intervals.^2)).^(1/2);
    %figure;
    %h = plot( fitresult4, x, y );
else %
    vq = Coeffs5(1)*intervals./(abs(Coeffs5(2)^2-intervals.^2)).^(1/2);
    %figure;
    %h = plot( fitresult5, x, y );
end

end

function [ vq, h ] = fitLeft( x, y, intervals)

[xData, yData] = prepareCurveData( x, y );

%Fit 1
% Set up fittype and options.
ft = fittype( 'abs(a/(x-b)^m)', 'independent', 'x', 'dependent', 'y' );
opts1 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts1.DiffMinChange = 1e-12;
opts1.Display = 'Off';
opts1.Lower = [-Inf min(x)*0.95 0.000001];
opts1.MaxIter = 1000;
opts1.Robust = 'Bisquare';
opts1.StartPoint = [rand min(x)*0.999999, 1];
opts1.TolFun = 1e-08;
opts1.Upper = [Inf min(x1)*0.999999999999 20];
% Fit model to data.
[fitresult1, gof] = fit( xData, yData, ft, opts1 );
Coeffs1 = coeffvalues(fitresult1);
rval1 = gof.rsquare;

%Fit 2
% Set up fittype and options.
ft = fittype( 'poly6' );
% Fit model to data.
[fitresult2, gof] = fit( xData, yData, ft );
Coeffs2 = coeffvalues(fitresult2);
rval2 = gof.rsquare;

%Fit 3
% Set up fittype and options.
ft = fittype( 'a*x/(abs(x^2*c-b^2))^(1/2)', 'independent', 'x', 'dependent', 'y' );
opts3 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts3.Algorithm = 'Levenberg-Marquardt';
opts3.Display = 'Off';
opts3.Robust = 'LAR';
inc=0;
rval3=0;
check=1;
while ( (rval3 < 0.96 || check>=0.05) && inc < 5 )
    opts3.StartPoint = [rand rand rand];
    % Fit model to data.
    [fitresult3, gof] = fit( xData, yData, ft, opts3 );
    rval3 = gof.rsquare;
    inc= inc+1;
    b = fitresult3.b;
    c = fitresult3.c;
    check = min(x)^2*c-b^2;
end
Coeffs3 = coeffvalues(fitresult3);

%Fit 4
% Set up fittype and options.
ft = fittype( 'a*x/(abs((x)^2*c-b^2))^(1/2)', 'independent', 'x', 'dependent', 'y' );
opts4 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts4.Display = 'Off';
opts4.Lower = [0 0 1];
opts4.Robust = 'LAR';
rval4=0;
inc=0;
while (rval4<0.93 && inc<5)
    opts4.StartPoint = [rand rand rand];
    opts4.Upper = [Inf min(x)*0.99999999 Inf];
    % Fit model to data.
    [fitresult4, gof] = fit( xData, yData, ft, opts4 );
    c = fitresult4.c;
    b = fitresult4.b;
    test = c*min(x)^2-b^2;
    if (test<0)
        rval4=0;
    else
        rval4 = gof.rsquare;
    end
    inc=inc+1;
end
Coeffs4 = coeffvalues(fitresult4);

%Fit 5
% Set up fittype and options.
str=strcat('a*x./(abs(x^2-',num2str(min(x)*0.9999999,12),'.^2))^(1/2)');
ft = fittype( str,'independent', 'x', 'dependent', 'y' );
opts5 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts5.Display = 'Off';
opts5.Lower = 0 ;
opts5.Robust = 'LAR';
opts5.Upper = Inf ;
rval5=0;
inc=0;
while (rval5<0.93 && inc<5)
    opts5.StartPoint = rand;
    % Fit model to data.
    [fitresult5, gof] = fit( xData, yData, ft, opts5 );
    rval5 = gof.rsquare;
    inc=inc+1;
end
Coeffs5 = coeffvalues(fitresult5);

if (rval3>=rval2 && rval3>=rval1 && rval3>=rval4 && rval3>=rval5)
    vq = Coeffs3(1)*intervals./(abs((min(x)*0.99999999)^2-intervals.^2)).^(1/2);
    %figure;
    %h = plot( fitresult3, x, y );
elseif (rval2>=rval3 && rval2>=rval1 && rval2>=rval4 && rval2>=rval5)
    vq = Coeffs2(7)+Coeffs2(6)*intervals+Coeffs2(5)*intervals.^2+Coeffs2(4)*intervals.^3+Coeffs2(3)*intervals.^4+Coeffs2(2)*intervals.^5+Coeffs2(1)*intervals.^6;
    %figure;
    %h = plot( fitresult2, x,y );
elseif (rval1>=rval3 && rval1>=rval2 && rval1>=rval4 && rval1>=rval5)
    vq = abs(Coeffs1(1)./(intervals-Coeffs1(2)).^Coeffs1(3));
    %figure;
    %h = plot( fitresult1, x,y);
elseif (rval4>=rval3 && rval4>=rval2 && rval4>=rval1 && rval4>=rval5)
    vq = Coeffs4(1)*intervals./(abs((min(Emid)*0.99999999)^2-intervals.^2)).^(1/2);
    %figure;
    %h = plot( fitresult4, x,y );
else
    vq = Coeffs5(1)*intervals./(abs((min(Emid)*0.9999999)^2-intervals.^2)).^(1/2);
    %figure;
    %h = plot( fitresult5, x,y );
end
end

function printDOS( n, m, En, DOS)

if (exist('DOS_Library','dir')==0)
    mkdir('DOS_Library');
end

filename = strcat('DOS_Library/DOS_',num2str(n),'_',num2str(m),'.txt');
if (exist(filename,'file')==0)
    fid = fopen(filename,'w');
    fprintf(fid,'Energy [ eV ] DOS [ 1/(eV Ang)]\n');
    for i=1:length(En)
        fprintf(fid,'%f \t %e\n',En(i),DOS(i));
    end
    fclose(fid);
end
end

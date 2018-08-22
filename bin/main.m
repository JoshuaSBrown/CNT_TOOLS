function main( n_values, m_values, Temperature_values, Resolution)

% This file is designed to calculate the carrier concentration of a range
% of CNTs

if (exist('CarrierConc','dir')~=7)
   mkdir('CarrierConc'); 
end

for index=1:length(n_values)
    n = n_values(index);

    for index2=1:length(m_values)
        m = m_values(index2);
        
        for index3=1:length(Temperature_values)
            Temperature = Temperature_values(index3);
            paste=0;
            num=0;
            if (exist('CarrierConc/CarrierConc.txt','file')==2)
                fid = fopen('CarrierConc.txt','r');
                [A] = fscanf(fid,'%d , %d %g [1/cm3] %g [1/cm3] %g [Ang]');
                [row, ~]=size(A);
                num = row/5;
                fclose(fid);
            end
            [carrierConc, carrierConc2, CNT_D] = DispBandDOS( n, m, Temperature, Resolution);
            
            fid = fopen('CarrierConc/CarrierConc.txt','w');
            if num==0
                fprintf(fid,'%d, %d %e [K] %e  [1/cm3]   %e  [1/cm3]   %f  [Ang]\n',n,m,Temperature,carrierConc, carrierConc2, CNT_D);
            else
                for i=1:num
                    
                    nn = A((i-1)*5+1);
                    mm = A((i-1)*5+2);
                    carrier = A((i-1)*5+3);
                    carrier2 = A((i-1)*5+4);
                    Diam = A((i-1)*5+5);
                    if ((nn==n && mm==m || nn>n && mm>m) && paste==0)
                        fprintf(fid,'%d, %d %e [K] %e  [1/cm3]   %e  [1/cm3]   %f  [Ang]\n',n,m,Temperature,carrierConc, carrierConc2, CNT_D);
                        paste=1;
                    end
                    fprintf(fid,'%d, %d %e [K] %e  [1/cm3]   %e  [1/cm3]   %f  [Ang]\n',nn,mm,Temperature,carrier, carrier2, Diam);
                end
            end
        end
        fclose(fid);
    end
    
end

end
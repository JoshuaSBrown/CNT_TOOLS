function main( chiral_vectors, Temperature_values, Resolution)

HelpMessage();

if (exist('CarrierConc','dir')~=7)
   mkdir('CarrierConc'); 
end

[rows, ~] = size(chiral_vectors);

for index=1:rows
    n = chiral_vectors(index,1);

    for index2=1:rows
        m = chiral_vectors(index2,2);
        
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

    function HelpMessage()
        fprintf('This script excepts the following flags followed by ');
        fprintf('thier associated values:\n\n');
        fprintf('--help [none]\n');
        fprintf('--chiral_vector [matrix(:,2)]\n');
        fprintf('--temperature [matrix(:,1)]\n');
        fprintf('--resolution [positive integer]\n\n');
        fprintf('--help displays acceptable flags and their appropraite ');
        fprintf('values\n\n');
        fprintf('--chiral_vector accepts a matrix with an arbitrary ');
        fprintf('number of rows and two columns. Columne 1 is associated');
        fprintf(' with the n chiral number of a CNT and column 2 is ');
        fprintf('associated with the m chiral number of a CNT.\n\n');
        fprintf('--temperature should be passed an array/matrix with ');
        fprintf('with a single column containing all the temperatures in');
        fprintf('Kelvin to be looped through.\n\n');
        fprintf('--resolution this should be passed a positive integer ');
        fprintf('value determing the resolution of the numerical ');
        fprintf('calculations used internally.');
    end
end
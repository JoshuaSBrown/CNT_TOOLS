function [ p, q ] = CNT_Symmetry_Vec_pq( n, m )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

%Initialize variables
q=-1;
p=1;
val=-1;
count=0;

%Calculations
[ t1, t2 ] = CNT_Translational_Vec_t1t2( n, m);

[ N ] = CNT_UnitCell_Num_Hex( n, m);

while((m*p-n*q<=0 || m*p-n*q>N || t1*q-t2*p~=1 || gcd(p,q)~=1) || count==0);
    count=1;
    q=-val;
    p=val;
  while (m*p-n*q<=-m*val-n*val && q >= val)
    
   
    while(m*p-n*q<=-m*val-n*val && p<=-val )
        
        %Check all Conditions
        if(t1*q-t2*p==1 && gcd(p,q)==1 && m*p-n*q>0 && m*p-n*q<=N)
            break
        end
        p=p+1;
    end
    
    if(t1*q-t2*p==1 && gcd(p,q)==1 && m*p-n*q>0 && m*p-n*q<=N)
        break
    end
    q=q-1;
    p=val;
    
  end
  
  if(t1*q-t2*p==1 && gcd(p,q)==1 && m*p-n*q>0 && m*p-n*q<=N)
    break
  end
  val=val-1;
end

%Post Check
assert(gcd(p,q)==1,'Greatest common devisor of p and q is not 1');
assert(t1*q-t2*p==1,'t1*q-t2*p does not equal 1');
assert(0<m*p-n*q,'m*p-n*q is not less than 0');
assert(N>=m*p-n*q,'m*p-n*q is not greater than N');

end


function [ L,U ] = LU_naive( A )

[m,n] = size(A);

if( m~= n)
    disp('error en las dimensiones de la matriz')
    return
end

U = A;
L = zeros(n,n);

for k = 1:1: n-1
    
    L(k,k) = 1.0;
        
    for i = k+1:1:n
        
        alfa = U(i,k)/ U(k,k);
        
        L(i,k) = alfa;
        U(i,k) = 0.0;
        
        for j = k+1:1:n % Fi = Fi - alfa*Fk
                                
            U(i,j) = U(i,j) - alfa*U(k,j);
                        
%           fprintf('k=%d i=%d j=%d alfa=%5.3f A(i,j)=%5.3f \n', k,i,j, alfa, A(i,j))
%           pause
            
        end
               
    end
        
end

L(n,n) = 1.0;

fprintf('error en la factorizacion LU = %f \n', norm(L*U-A, inf))

%U = triu(A);
%L = tril(A,-1) + eye(n);

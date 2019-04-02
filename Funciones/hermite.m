function [s, z] = hermite(x, y, y1)
    [n, m] = size(x);
    s = NaN*ones(2*n,1);
    q = NaN*zeros(2*n, 2*n);
    z = NaN*zeros(2*n,1);
    z(1:2:2*n) = x;
    z(2:2:2*n) = x;
    k = 1;
    for i = 1:2*n
        q(i, 1) = y(k);
        if mod(i,2)==0
            k = k+1;
        end
    end
    s(1) = y(1);
    for i = 3:2:2*n-1
        q(i,2) = (q(i, 1)-q(i-1, 1))/(z(i)-z(i-1));
    end
    s(2) = y1(1);
    for i = 2:2:2*n
        q(i,2) = y1(i/2);
    end
    for j = 3: 2*n
        for i = j:2*n
            q(i,j) = (q(i, j-1)-q(i-1, j-1))/(z(i)-z(i-j+1));
        end
        s(j) = q(j,j);
    end
        
end
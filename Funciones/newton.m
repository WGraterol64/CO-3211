function s = newton(x, y)
  [n,m] = size(x);
  s = NaN*ones(n,1);
  q = zeros(n, n);
  for i = 1:n
    q(i, 1) = y(i);
  end
  s(1) = y(1);
  for j = 2: n
    for i = j:n
      q(i,j) = (q(i, j-1)-q(i-1, j-1))/(x(i)-x(i-j+1));
    end
    s(j) = q(j,j);
  end
endfunction

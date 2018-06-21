function x = rayleigh(A, epsilon, mu, x)
  x = x / norm(x);
  % the backslash operator in Octave solves a linear system
  y = (A - mu * eye(A(:).')) \ x; 
  lambda = y' * x;
  mu = mu + 1 / lambda
  err = norm(y - lambda * x) / norm(y)

  while err > epsilon
    x = y / norm(y);
    y = (A - mu * eye(A(:).')) \ x;
    lambda = y' * x;
    mu = mu + 1 / lambda
    err = norm(y - lambda * x) / norm(y)
  end

end
function y = HG(n, q, k, x )
%HG gives the Hermite-Gauss function of order n.
%   One-dimensional normalized Hermite-Gauss mode of order n with and
%   complex beam parameter q.
%   
%   Inputs arguments:
%   n 	- Integer. Order of Hermite polynomial.
%   q 	- Complex Scalar. Complex beam parameter = z-izR.
%	k 	- Scalar. Wavenumber.
%   x 	- Array. Function argument.


z = real(q);
w = abs(sqrt(2/imag(k/q)));


y = mfun('H', n, sqrt(2).* x ./ w) .* exp(1i*k*x.^2/(2*q) + 1i*k*z );
y = y ./ abs(sqrt(sum(y.^2)));

end


function [rho, omega2] = Gbar_Eliashberg(lambdaV, Ncut, mustar)
	%% Determine omega2 by extrapolating n^2 * lambda(n)
	nw = Ncut + 1;
	ns = (nw-3) : nw;
	w2fit = fit(1./(ns'-1).^2, ((ns-1).^2 .* lambdaV(ns))', 'poly1');
	omega2 = sqrt(w2fit.p2 / lambdaV(1)) * 2 * pi;  % in unit of kB*T

	mustarN = 1/(1/mustar + log(omega2/(pi*(2*Ncut+1))));

	K = zeros(Ncut+1, Ncut+1);

	for m = 0:Ncut
		for n = 0:m
			K(m+1, n+1) = lambdaV(m-n+1) - 2*mustarN;
			if m+n+2 <= nw
				K(m+1, n+1) = K(m+1, n+1) + lambdaV(m+n+2);
			else
				K(m+1, n+1) = K(m+1, n+1) + lambdaV(1) * omega2^2 / (2*pi*(m+n+1))^2;
			end

			if m == n
	  			K(m+1, m+1) = K(m+1,m+1) - (2*m+1) - lambdaV(1) ...
	  						- 2 * sum(lambdaV(2:(m+1)));
	  		end

			K(n+1, m+1) = K(m+1, n+1);
		end
	end

	eigK = eig(K);

	rho = max(eigK);
end
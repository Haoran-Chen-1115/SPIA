function [chitau1, omegatau1, hqtau1] = extrapolatePDF(NBEAD, NBEAD1, omegas, hqs)

	nq = size(omegas, 1);
	omegatau = real(fft(omegas.'));
	omegatau1 = zeros(nq, NBEAD1);
	tau = (0:(NBEAD-1))/NBEAD;
	tau1 = (0:(NBEAD1-1))/NBEAD1;
	for k = 1:nq
		id = omegatau(:, k) > 0;
		logomegatau = @(t) interp1([tau(id), 1], [log(omegatau(id, k)); log(omegatau(1, k))], ...
				t, 'pchip');
		omegatau1(k, :) = exp(logomegatau(tau1)).';
	end

	hqtau = real(fft(hqs.'));
	hqtauf = @(t) interp1([tau, 1], [hqtau; hqtau(1, :)], t);
	hqtau1 = hqtauf(tau1).';

	chitau1 = omegatau1 + hqtau1;
end

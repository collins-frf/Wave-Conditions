function [freq, direction, spec1] = TMA_2peak(fmin,fmax,df,FM1,gam1,Hs1,nn1,theta1,FM2,gam2,Hs2,nn2,theta2)
%
%fmin = minimum spectral freuency
%fmax = maximum spectral frequency
%df = freguency resolution

% first spectral peak (second spectral peak is *2)
%FM1 = peak frequency
%gam1 = frequency spreading
%Hs1 = significant wave height
%nn1 = directional spreading
%theta1 = peak direction


hm1 = 500; hm2=hm1;

freq = [fmin:df:fmax];
NF = length(freq);
da = 5; % hardwired to 5 deg resolution in STWAVE
direction = [0:da:(theta1+90)]'; % hardwired in STWAVE
NA = length(direction); 

%calculate the TMA spectrum
f = repmat(freq,NA,1);
a = repmat(direction,1,NF);
omh1 = 2*pi*sqrt(hm1/9.8)*f;
if omh1<1
  phi1 = 0.5*(omh1.^2);
elseif omh1>=1 & omh1<=2
  phi1 = 1 - 0.5*((2 - omh1).^2);
else
  phi1 = 1;
end
sigma = 0.07*ones(size(f));
sind1 = find(f>=FM1);
sigma(sind1) = 0.09;
specF1 = (9.8^2)./((f.^5)*((2*pi).^4)).*phi1.*exp((-1.25*(FM1./f).^4) + log(gam1)*exp(-((f - FM1).^2)./(2*(sigma.^2)*(FM1^2))));
specA1 = cos((a-theta1)*pi/180).^nn1;
spec1 = specA1.*specF1; 
spec1(isnan(spec1))=0;
normSpec1 = sum(spec1(:))*df*da/(Hs1^2/16);
spec1 = spec1/normSpec1;
% spec.in  -SECOND
%calculate the TMA spectrum
omh2 = 2*pi*sqrt(hm2/9.8)*f;
if omh2<1
  phi2 = 0.5*(omh2.^2);
elseif omh2>=1 & omh2<=2
  phi2 = 1 - 0.5*((2 - omh2).^2);
else
  phi2 = 1;
end
sind2 = find(f>=FM2);
sigma(sind2) = 0.09;
specF2 = (9.8^2)./((f.^5)*((2*pi).^4)).*phi2.*exp((-1.25*(FM2./f).^4) + log(gam2)*exp(-((f - FM2).^2)./(2*(sigma.^2)*(FM2^2))));
specA2 = cos((a-theta2)*pi/180).^nn2;
spec2 = specA2.*specF2;
spec2(isnan(spec2))=0;
normSpec2 = sum(spec2(:))*df*da/(Hs2^2/16);
spec2 = spec2/normSpec2;
spec=(spec1+spec2);

save("TMA_out")

subplot(1,3,1)
sgtitle("Hs = " + Hs1 + "m, Dir = " + theta1 + "*, F = " + FM1 + "Hz");
s = pcolor(freq, direction, spec1);
scb = colorbar;
title('2D Wave Spectra');
xlabel('Frequency (Hz)');
ylabel('Direction (0 is E)');
s.EdgeColor = 'none';
s.FaceColor = 'interp';

spec_freq_mean = mean(spec1,1);
subplot(1,3,2)
s = plot(freq, spec_freq_mean);
title('Avg. Wave Energy at each Frequency');
xlabel('Frequency (Hz)');
ylabel('Avg. Wave Energy');

spec_dir_mean = mean(spec1,2);
subplot(1,3,3)
s = plot(direction, spec_dir_mean);
title('Avg. Wave Energy at each Direction');
xlabel('Direction (*)');
ylabel('Avg. Wave Energy');





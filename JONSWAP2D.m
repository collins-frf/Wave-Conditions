% Generate a JONSWAP2D spectrum (freq and direction) from the given wave parameters
% fmin = minimum spectral frequency
% fmax = max spectral frequency
% df = frequency resolution
% dtheta = directional resolution
% fpeak = peak spectral frequency
% thp = peak direction
% Hs = significant wave height
% gam = frequency spread
% nn = directional spread

function [f,SfG,theta,G]=JONSWAP2D(fmin,fmax,df,dtheta,fpeak,thp,Hs,gam,nn)

f=fmin:df:fmax;

Tp=1/fpeak; fp = fpeak;
Bj=0.0624*(1.094-0.01915*log(gam))/(.230+0.0336*gam-0.185*((1.9+gam)^-1));

for j=1:length(f);
if f(j) < fpeak;
sigma(j)=0.07;
else
sigma(j)=0.09;
end
end

Sf=Bj.*(Hs.^2).*(Tp.^-4).*(f.^-5).*exp(-1.25.*((Tp.*f).^-4)).*gam.^(exp(-((Tp.*f-1).^2)./(2.*sigma.^2)));

smax=600;
theta=[-90:dtheta:90].*pi/180;
tp = thp*pi/180;
%for j=1:length(f);
%if f(j) < fp;
%s(j)=((f(j)/fp)^5)*smax;
%else
%s(j)=((f(j)/fp)^(-2.5))*smax;
%end
%end


%G0=(1./pi).*(2.^(2.*s-1)).*(((gamma(s+1)).^2)./(gamma(2.*s+1)));

for i=1:length(theta)
G(i,:)=cos(theta(i)-tp).^(nn);
end

for kk=1:length(theta)
SfG(kk,:)=fliplr(Sf(1,:)).*G(kk,:);
end

figure
surf(fliplr(f),theta,SfG); shading interp

dspr = sqrt(sum(G'.*(2*sin(theta./2)).^2)*dtheta)*180/pi

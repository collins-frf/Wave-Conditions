function [Hs,emean,Tp,Tm,f,Se,f_ave,S_ave]=spectral(t,e)
N=length(t);
delt=mean(diff(t));
emean=mean(e);
esub = e-emean;

% compute Fourier coefficients for eta
Ce = fft(esub,N);
% estimate spectrum for u
Se = Ce.*conj(Ce)/N;
Se(round(N/2)+1+1:N)= [ ];
delf = 1/(delt*N);
Se(2:round(N/2)) = 2*Se(2:round(N/2))/(N*delf);
% check
delf;
var3=sum(Se)*delf;
% calculate Hm0
Hs = 4.004*sqrt(var3);
f=delf*[0:length(Se)-1];

f_ave=0;
S_ave=0;

maxSe=max(Se);
Tp=0;
for i=1:length(Se)
   if Se(i)==maxSe
      Tp=1/f(i);
   end
end

var4=sum(Se.*f)*delf;
if Tp==0
   Tm=0;
else
   Tm=var3/var4;
end

%save spectrum_out.mat



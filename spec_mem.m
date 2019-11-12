function D = spec_mem(r1,r2,alpha1,alpha2,dirs)

%MEM Maximum Entropy Method calculation for directional distribution function
%   D = mem(R1,R2,ALPHA1,ALPHA2,DIRS) calculates the directional distribution
%   function for a 2D spectrum for NDBC parameters using the Maximum Entropy
%   Method as outlined by Lygre & Krogstad (1986 - JPO). The inputs are 
%   defined on http://www.ndbc.noaa.gov/measdes.shtml - R1 and R2 are the 
%   two spreading parameters (NOTE: 0<r1,r2<1), ALPHA1 and ALPHA2 are 
%   the mean and principal wave directions, and DIRS is the directional 
%   space of the distribution.
%
%   D = mem(r1,r2,alpha1,alpha2,dirs)
%
% Dave Thompson (dthompson@usgs.gov)

%% Setup to handle f and t!
deg1 = ((270-alpha1)<0)*360 + (270-alpha1);
deg2 = ((270-alpha2)<0)*360 + (270-alpha2);
dirs = ((270-dirs)<0)*360 + (270-dirs);

% Calculate the fourier coefficients.
a1 = r1.*cosd(deg1);
b1 = r1.*sind(deg1);
a2 = r2.*cosd(2*deg2);
b2 = r2.*sind(2*deg2);

c1 = a1 + 1i*b1;
c2 = a2 + 1i*b2;

phi1 = (c1-c2.*conj(c1))./(1-abs(c1).^2);
phi2 = c2-c1.*phi1;

num = real(1 - phi1.*conj(c1) - phi2.*conj(c2));
deg = dirs*pi/180;

% assuming r1, r2, alpha1, alpha2 came in as [t,f];
% phi1, phi2 and num are now [t,f]. Need [f,d,t];
phi1 = permute(repmat(phi1.',[1 1 length(deg)]),[1 3 2]);
phi2 = permute(repmat(phi2.',[1 1 length(deg)]),[1 3 2]);

%phi1 = repmat(phi1.',[1 length(deg)]);
%phi2 = repmat(phi2.',[1 length(deg)]);
deg = repmat(deg,[size(phi1,1) 1 size(phi1,3)]);
den = abs((1 - phi1.*exp(-1i*deg) - phi2.*exp(-1i*2*deg)).^2);

num = permute(repmat(num',[1 1 size(deg,2)]),[1 3 2]);
%num = repmat(num',[1 size(deg,2)]);
D = num./(2*pi*den);
D = D./repmat(sum(D,2),[1 size(D,2) 1]);
%D = D/sum(D);
return
%% Original setup (per f)
% Transform directions from Nautical (CW from N) to Cartesian (CCW from E).
deg1 = ((270-alpha1)<0)*360 + (270-alpha1);
deg2 = ((270-alpha2)<0)*360 + (270-alpha2);
dirs = ((270-dirs)<0)*360 + (270-dirs);

% Calculate the fourier coefficients.
a1 = r1*cosd(deg1);
b1 = r1*sind(deg1);
a2 = r2*cosd(2*deg2);
b2 = r2*sind(2*deg2);

%%
% Maximum Entropy Method - Lygre & Krogstad (1986 - JPO)
% Eqn. 13:
% phi1 = (c1 - c2c1*)/(1 - abs(c1)^2)
% phi2 = c2 - c1phi1
% 2piD = (1 - phi1c1* - phi2c2*)/abs(1 - phi1exp(-itheta) -phi2exp(2itheta))^2
%     
% phi1 & phi2 come from the Yule-Walker equations (Eqn. 5):
% -       - -      -   -    -
% | 1 c1* | | phi1 | = | c1 |  gives: phi1 + phi2c1* = c1
% | c1 1  | | phi2 |   | c2 |         phi1c1 + phi2 = c2
% -       - -      -   -    - 
%
% c1 and c2 are the complex fourier coefficients
c1 = a1 + 1i*b1;
c2 = a2 + 1i*b2;

phi1 = (c1-c2*conj(c1))/(1-abs(c1)^2);
phi2 = c2-c1*phi1;

num = real(1 - phi1*conj(c1) - phi2*conj(c2));
deg = dirs*pi/180;
den = abs((1 - phi1*exp(-1i*deg) - phi2*exp(-1i*2*deg)).^2);

D = num./(2*pi*den);
D = D/sum(D);
return

%% Setup to handle all f's
deg1 = ((270-alpha1)<0)*360 + (270-alpha1);
deg2 = ((270-alpha2)<0)*360 + (270-alpha2);
dirs = ((270-dirs)<0)*360 + (270-dirs);

% Calculate the fourier coefficients.
a1 = r1.*cosd(deg1);
b1 = r1.*sind(deg1);
a2 = r2.*cosd(2*deg2);
b2 = r2.*sind(2*deg2);

c1 = a1 + 1i*b1;
c2 = a2 + 1i*b2;

phi1 = (c1-c2.*conj(c1))./(1-abs(c1).^2);
phi2 = c2-c1.*phi1;

num = real(1 - phi1.*conj(c1) - phi2.*conj(c2));
deg = dirs*pi/180;

phi1 = repmat(phi1.',[1 length(deg)]);
phi2 = repmat(phi2.',[1 length(deg)]);
deg = repmat(deg,[size(phi1,1) 1]);
den = abs((1 - phi1.*exp(-1i*deg) - phi2.*exp(-1i*2*deg)).^2);

num = repmat(num',[1 size(deg,2)]);
D = num./(2*pi*den);
D = D./repmat(sum(D,2),[1 size(D,2)]);
return



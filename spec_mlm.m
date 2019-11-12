function D = spec_mlm(r1,r2,alpha1,alpha2,dirs)

%MLM Maximum Likelihood Method calculation for directional distribution function
%   D = mlm(R1,R2,ALPHA1,ALPHA2,DIRS) calculates the directional distribution
%   function for a 2D spectrum for NDBC parameters using the Maximum
%   Likelihood Method as outlined by Oltman-Shay & Guza (1984 - JPO). The 
%   inputs are defined on http://www.ndbc.noaa.gov/measdes.shtml - R1 and 
%   R2 are the two spreading parameters (NOTE: 0<r1,r2<1), ALPHA1 and 
%   ALPHA2 are the mean and principal wave directions, and DIRS is the 
%   directional space of the distribution.
%
%   D = mlm(r1,r2,alpha1,alpha2,dirs)
%
% Dave Thompson (dthompson@usgs.gov)

%%
% Set k and c11 to unity. 
k = 1;
c11 = 1;

% Transform directions from Nautical (CW from N) to Cartesian (CCW from E).
dirs = ((270-dirs)<0)*360 + (270-dirs);

% Calculate co and quad-spectrum (from Earle, 1999)
c22 = (1/2)*k^2*c11*(1-r2*cosd(2*alpha2));
c33 = (1/2)*k^2*c11*(1+r2*cosd(2*alpha2));
c23 = (1/2)*k^2*c11*r2*sind(2*alpha2);
q12 = -k*c11*r1*sind(alpha1);
q13 = -k*c11*r1*cosd(alpha1);
c32 = conj(c23);
q21 = -q12;
q31 = -q13;

% cross-spectral data matrix
M = [c11 1i*q12 1i*q13;
   1i*q21 c22 c23;
   1i*q31 c32 c33];
M = inv(M); % for Eq. 6 in O-S & G.

% G
for ii=1:length(dirs)
   G(ii,1) = 1;
   G(ii,2) = 1i*cosd(dirs(ii));
   G(ii,3) = 1i*sind(dirs(ii));
end

% Calculate direction distribution -> Eq. 6 in O-S & G.
D = zeros(1,length(dirs));
for ii=1:length(dirs)
   for n=1:3
      for m=1:3
         D(ii) = D(ii) + M(n,m) * G(ii,n) * conj(G(ii,m));
      end
   end
   D(ii) = 1/D(ii);
end
D = D/sum(D);
return
%% ORIGINAL 
% I got this fortran program from Bill
% This program is to use the MLM to compute the wave directional
% distribution function based on the given r1,r2,alph1, and alph2.
%
% firt convert the r1,r2,alp1,alp2, to 
% obtain c22,c33,c12,c13,q12,q13. (c11 was set up as one)
%
% because in this program the co and quad are used to obtain direction 
% distribution fuction for a given frequency the direction distributioin 
% will be normalized to unity later so we can set the xk and c11 to be 
% unity and will not affect the end result. because this unity set up for 
% c11 and xk, the c22, c33, c23, q12 and q13 obtained here are only good 
% for this calculation should not be considered for other purposes.

xk=1; % wave number

% ref numbers are equation numbers in M.D. Earle et. al. (1999)
c11=1;                        % (29)
tmp1=(xk*xk)/2.;              % (30)
tmp6=cos(2*alpha2*pi/180);    
tmp7=sin(2*alpha2*pi/180);
tmp8=sin(alpha1*pi/180);
tmp9=cos(alpha1*pi/180);
c22=c11*tmp1*(1-r2*tmp6);     % (31)
c33=c11*tmp1*(1+r2*tmp6);     % (32)
c23=c11*tmp1*r2*tmp7;         % (33)
q12=-c11*xk*r1*tmp8;          % (34)
q13=-c11*xk*r1*tmp9;          % (35)
xk=((c22+c33)/c11).^0.5;      % (18)
d=c11*c22*c33-c11*(c23.^2)-c33*(q12^2)-c22*(q13^2)+2.*q12*q13*c23;
if(abs(d) >= 0.0000000001)
   temp=2*(c22*c33-c23^2)/d;
   a0=temp+(c11*c22+c11*c33-q12^2-q13^2)*(xk^2)/d;
   a1=-2*(q12*c33-q13*c23)*xk/d;
   b1=2.*(q12*c23-q13*c22)*xk/d;
   a2=0.5*(c11*c33-c11*c22-q13^2+q12^2)*(xk^2)/d;
   b2=-(c11*c23-q12*q13)*(xk^2)/d;
   %[a0 a1 a2 b1 b2]
   ndd=270-wavedir;
   denim=0.5*a0+a1*cos(ndd*pi/180)+b1*sin(ndd*pi/180)+...
      a2*cos(2*ndd*pi/180)+b2*sin(2*ndd*pi/180); % (6)
   dsprd=abs(denim.^(-1));
   sum1=sum(dsprd);
   dirdismlm=dsprd/sum1;
else
   dirdismlm=wavedir*0;
end
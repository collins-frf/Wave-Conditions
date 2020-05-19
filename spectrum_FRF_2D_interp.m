function spectrum_FRF_2D_interp(f,theta,E_D,Hs_o,Tp_o,n_cutoff)
% % % Interpolates the 2D frequency x direction spectrum 
% from FRF conventions of 62 by 75 to a value that is 
% most appropriate for model setup. 
% INPUT
%   f - input wave frequencies shaped [nf]
%   theta -- input wave directions shaped [nd]
%   E_D - input spectra shaped [nf, nd]
%   Hs_o - wave height 
%   Tp_o - wave peak period 
%   n_cutoff - the frequency cutoff by number of grid points in X and Y, Current limitation of the model
% written by Pat Lynette, modified by spicer bak
rng('shuffle');  % random rand seed
RunDurationMin = 10; % duration in minutes of simulation
min_period=3.5;  % min allowable period
min_theta=-90; % min allowable theta
max_theta=90; % max allowable theta

Hmo = 0;
nf=length(f);

for i=1:nf
   for j=1:length(theta)
      if i==1
         del_f=(f(2)-f(1));
      elseif i==length(f)
         del_f=(f(length(f))-f(length(f)-1));
      else
         del_f=(f(i+1)-f(i))/2.+(f(i)-f(i-1))/2.;
      end
      
      if j==1
         del_theta=(theta(2)-theta(1));
      elseif j==length(theta)
         del_theta=(theta(length(theta))-theta(length(theta)-1));
      else
         del_theta=(theta(j+1)-theta(j))/2.+(theta(j)-theta(j-1))/2.;
      end
      
      Hmo = Hmo + E_D(i,j)*del_f*del_theta;
   end
end

Hmo_full_spectrum = sqrt(Hmo)*4.004;
% find min and max values for truncation
f_max=1/min_period;
f_max_ind=find(f>f_max,1);
t_min_ind=find(theta>min_theta,1);
t_max_ind=find(theta>max_theta,1);
% actually truncate 
f=f(1:f_max_ind);
theta=theta(t_min_ind:t_max_ind);
E_D=E_D(1:f_max_ind,t_min_ind:t_max_ind);
% calculate number of values left
nf=length(f);
nt=length(theta);
%%
Hmo_c=0;
for i=1:nf
    for j=1:length(theta)
        if i==1
            del_f=(f(2)-f(1));
        elseif i==length(f)
            del_f=(f(length(f))-f(length(f)-1));
        else
            del_f=(f(i+1)-f(i))/2.+(f(i)-f(i-1))/2.;
        end
        
        if j==1
            if length(theta)==1
                del_theta=1;
            else
                del_theta=(theta(2)-theta(1));
            end
        elseif j==length(theta)
            del_theta=(theta(length(theta))-theta(length(theta)-1));
        else
            del_theta=(theta(j+1)-theta(j))/2.+(theta(j)-theta(j-1))/2.;
        end
        
        Hmo_c = Hmo_c + E_D(i,j)*del_f*del_theta;
    end
end
Hmo_cut_spectrum = sqrt(Hmo_c)*4.004;

% refine df
df_new=1/(RunDurationMin*60);  
f_new=[f(1):df_new:f_max];

repeat_time=1/df_new/60; % cycle repeat time in minutes
disp(['Time series cycle time (min): ' num2str(repeat_time)])

dt_new=16;
theta_new=[theta(1):dt_new:theta(nt)];
E_D=interp2(f,theta,E_D',f_new,theta_new','spline');
E_D=E_D';
f=f_new;
nf=length(f);
theta=theta_new;
nt=length(theta);
%%
E_Dsort=sort(reshape(E_D,[nf*nt,1]),'descend');
E_cutoff=E_Dsort(min(length(E_Dsort),n_cutoff));

Hmo_c2=0;
for i=1:nf
    for j=1:length(theta)
        if i==1
            del_f=(f(2)-f(1));
        elseif i==length(f)
            del_f=(f(length(f))-f(length(f)-1));
        else
            del_f=(f(i+1)-f(i))/2.+(f(i)-f(i-1))/2.;
        end
        
        if j==1
            del_theta=(theta(2)-theta(1));
        elseif j==length(theta)
            del_theta=(theta(length(theta))-theta(length(theta)-1));
        else
            del_theta=(theta(j+1)-theta(j))/2.+(theta(j)-theta(j-1))/2.;
        end
        
        if E_D(i,j)>=E_cutoff
            Hmo_c2 = Hmo_c2 + E_D(i,j)*del_f*del_theta;
        end
    end
end
Hmo_cut_spectrum2 = sqrt(Hmo_c2)*4.004;   % calculate the cut spectrum wave height
disp(['Input ','Full Spectrum','Cut Spectrum'])
disp([Hs_o, Hmo_full_spectrum, Hmo_cut_spectrum])

% adjust E_D
E_D=E_D*Hmo_c/Hmo_c2;  % should lead to top 900 E_D points yielding Hmo_cut_spectrum
E_cutoff=E_cutoff*Hmo_c/Hmo_c2;

fileID = fopen('irrWaves.txt','w');
line1='\n';
line2=['[NumberOfWaves] ' num2str(n_cutoff) '\n'];
line3='=================================\n';

for i=1:3
    eval(['cline = line' num2str(i) ';'])
    fprintf(fileID,cline);	
end

% adjust factor - model underpredicts generation of shorter waves
Tfact=[5 6];
adfac=[1.1 1.]*1.1;

cur_ind=0;
for i=1:nf
    if i==1
        del_f=(f(2)-f(1));
    elseif i==nf
        del_f=(f(nf)-f(nf-1));
    else
        del_f=(f(i+1)-f(i))/2.+(f(i)-f(i-1))/2.;
    end
    
    Tc=1/(f(i));
    if Tc<Tfact(1)
        adjust_c=adfac(1);
    elseif Tc>Tfact(2)
        adjust_c=adfac(2);
    else
        adjust_c=interp1(Tfact,adfac,Tc);
    end
    
    for j=1:length(theta)
        if j==1
            del_theta=(theta(2)-theta(1));
        elseif j==length(theta)
            del_theta=(theta(length(theta))-theta(length(theta)-1));
        else
            del_theta=(theta(j+1)-theta(j))/2.+(theta(j)-theta(j-1))/2.;
        end
        
        if E_D(i,j)>=E_cutoff
            cur_ind=cur_ind+1;
            if cur_ind>n_cutoff
                error(['Too many frequencies, something not correct in logic'])
            end
            amp(i,j)=sqrt(2.*E_D(i,j)*del_f*del_theta);
            file_data(cur_ind,1)=adjust_c*amp(i,j);
            file_data(cur_ind,2)=1/(f(i));
            file_data(cur_ind,3)=(270-(theta(j)+18.2))*3.1415/180;  % add 18.2 to account for site rotation
            file_data(cur_ind,4)=rand*2*3.1415;
            cline=[num2str(file_data(cur_ind,:)) '\n'];
            fprintf(fileID,cline);
        end
    end
end


fclose(fileID);

figure(3)
clf
subplot(2,1,2), plot(f,sum(E_D,2)*mean(diff(theta)));
ylabel('Energy (m^2s)','FontSize',5)
xlabel('Frequency (Hz)','FontSize',5)
title(['\Theta-Integrated Energy Spectrum'],'FontSize',5)
axis([-Inf Inf -Inf Inf])
set(gca,'fontsize',5)
subplot(2,1,1), pcolor(f,theta,E_D')
shading interp
title(['2-D Energy Spectrum, H_s (m) = ' num2str(Hmo_cut_spectrum) ', T_p (s )= ' num2str(Tp_o)],'FontSize',5)
ylabel(['Direction (degrees)'],'FontSize',5)
%xlabel(['Frequency (Hz)'],'FontSize',5)
axis([-Inf Inf 0 Inf])
set(gca,'fontsize',5)
%% output 
save repeat_time.txt repeat_time -ascii
% also make s input spectral file

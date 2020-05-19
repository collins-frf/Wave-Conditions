function spectrum_FRF(f,theta,E_D,Hs_o,Tp_o)

rng('shuffle');  % random rand seed

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


min_period=6;  % min allowable period
min_theta=0; % min allowable theta
max_theta=160; % max allowable theta

f_max=1/min_period;
f_max_ind=find(f>f_max,1);
t_min_ind=find(theta>min_theta,1);
t_max_ind=find(theta>max_theta,1);

f=f(1:f_max_ind);
theta=theta(t_min_ind:t_max_ind);
E_D=E_D(1:f_max_ind,t_min_ind:t_max_ind);

nf=length(f);

Hmo=0;
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
Hmo_cut_spectrum = sqrt(Hmo)*4.004;
disp(['Input ','Full Spectrum','Cut Spectrum'])
disp([Hs_o, Hmo_full_spectrum, Hmo_cut_spectrum])


fileID = fopen('irrWaves.txt','w');
line1='\n';
line2=['[NumberOfWaves] ' num2str(length(f)*length(theta)) '\n'];
line3='=================================\n';

for i=1:3
    eval(['cline = line' num2str(i) ';'])
    fprintf(fileID,cline);	
end

Hmo_truncated_spectrum=0;
cur_ind=0;
for i=1:nf
    if i==1
        del_f=(f(2)-f(1));
    elseif i==nf
        del_f=(f(nf)-f(nf-1));
    else
        del_f=(f(i+1)-f(i))/2.+(f(i)-f(i-1))/2.;
    end
    
    for j=1:length(theta)
        if j==1
            del_theta=(theta(2)-theta(1));
        elseif j==length(theta)
            del_theta=(theta(length(theta))-theta(length(theta)-1));
        else
            del_theta=(theta(j+1)-theta(j))/2.+(theta(j)-theta(j-1))/2.;
        end
        
        cur_ind=cur_ind+1;
        amp(i,j)=sqrt(2*E_D(i,j)*del_f*del_theta);
        file_data(cur_ind,1)=1.5*amp(i,j);
        file_data(cur_ind,2)=1/(f(i));
        file_data(cur_ind,3)=(270-(theta(j)+18.2))*3.1415/180;  % add 18.2 to account for site rotation
        file_data(cur_ind,4)=rand*2*3.1415;
        cline=[num2str(file_data(cur_ind,:)) '\n'];
        fprintf(fileID,cline);
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
subplot(2,1,1), pcolor(f,theta,amp')
shading interp
title(['2-D Amplitude Spectrum, H_s (m) = ' num2str(Hmo_cut_spectrum) ', T_p (s )= ' num2str(Tp_o)],'FontSize',5)
ylabel(['Direction (degrees)'],'FontSize',5)
%xlabel(['Frequency (Hz)'],'FontSize',5)
axis([-Inf Inf 0 Inf])
set(gca,'fontsize',5)


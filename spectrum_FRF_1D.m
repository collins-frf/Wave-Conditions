function spectrum_FRF_1D(f,theta,E_D,Hs_o,Tp_o)

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
nt=length(theta);
E_1D=zeros(nf,1);
D_1D=E_1D;

for i=1:nf
    ED_sum=0;
    for j=1:length(theta)
        
        if j==1
            del_theta=(theta(2)-theta(1));
        elseif j==length(theta)
            del_theta=(theta(length(theta))-theta(length(theta)-1));
        else
            del_theta=(theta(j+1)-theta(j))/2.+(theta(j)-theta(j-1))/2.;
        end
        
        E_1D(i)=E_1D(i)+ E_D(i,j)*del_theta;
        ED_sum=ED_sum+E_D(i,j)*theta(j);
    end
    D_1D(i)=ED_sum/sum(E_D(i,:));
end

% refine df
nf_new=750;
df_new=(f_max-f(1))/(nf_new-1);
1/df_new/60
f_new=[f(1):df_new:f_max];
E_D=interp1(f,E_1D,f_new,'spline');
theta_D=interp1(f,D_1D,f_new,'spline');
f=f_new;
nf=length(f);

Hmo=0;
for i=1:nf
    if i==1
        del_f=(f(2)-f(1));
    elseif i==length(f)
        del_f=(f(length(f))-f(length(f)-1));
    else
        del_f=(f(i+1)-f(i))/2.+(f(i)-f(i-1))/2.;
    end
    
    Hmo = Hmo + E_D(i)*del_f;
end
Hmo_cut_spectrum = sqrt(Hmo)*4.004;
disp(['Input ','Full Spectrum','Cut Spectrum'])
disp([Hs_o, Hmo_full_spectrum, Hmo_cut_spectrum])


fileID = fopen('irrWaves.txt','w');
line1='\n';
line2=['[NumberOfWaves] ' num2str(length(f)) '\n'];
line3='=================================\n';

for i=1:3
    eval(['cline = line' num2str(i) ';'])
    fprintf(fileID,cline);	
end

Hmo_truncated_spectrum=0;
for i=1:length(f)
    if i==1
        del_f=(f(2)-f(1));
    elseif i==length(f)
        del_f=(f(length(f))-f(length(f)-1));
    else
        del_f=(f(i+1)-f(i))/2.+(f(i)-f(i-1))/2.;
    end
    amp(i)=sqrt(2*E_D(i)*del_f);
    file_data(i,1)=amp(i);
    file_data(i,2)=1/f(i);
    file_data(i,3)=(270-(theta_D(i)+18.2))*3.1415/180;
    file_data(i,4)=rand*2*3.1415;
    cline=[num2str(file_data(i,:)) '\n'];
    fprintf(fileID,cline);
end

fclose(fileID);

figure(3)
clf
subplot(2,1,1), plot(1./f,E_D);
ylabel('Energy (m^2s)','FontSize',5)
xlabel('Wave Period (seconds)','FontSize',5)
title(['Energy Spectrum from CDIP'],'FontSize',5)
axis([-Inf Inf 0 Inf])
set(gca,'fontsize',5)
subplot(2,1,2), plot(1./f,theta_D);
ylabel('Mean Wave Direction (degrees)','FontSize',5)
xlabel('Wave Period (seconds)','FontSize',5)
title(['Direction from CDIP'],'FontSize',5)
set(gca,'fontsize',5)
axis([-Inf Inf -Inf Inf])





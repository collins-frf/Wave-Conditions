function [] = binBuoyData(waves)
%% function [] = binBuoyData(waves)
%%
%% INPUTS:  
%%	  waves - a structure generated by calling ndbc_get.m
%%
%%
%% Dave Thompson, Joe Long, Rangley Mickey - USGS 


B.time = waves.stdmet.time;
B.H = waves.stdmet.wave_height;
B.T = waves.stdmet.dominant_wpd;
B.D = waves.stdmet.mean_wave_dir;
B.ws = waves.stdmet.wind_spd;
B.wd = waves.stdmet.wind_dir;

fnames = fieldnames(B);

% find non NaN values/times
id = find(~isnan(B(1).H) & ~isnan(B(1).T) & ~isnan(B(1).D) & B(1).D<999);

for jj=1:length(fnames)
   B.(fnames{jj}) = B.(fnames{jj})(id);
end
N = length(id);

% right now hard-code the bin sizes...new places in the world (e.g. not the Gulf of Mexico) may change this
Hbins = [0:.5:2 3 4 5 6];
Hshift = 0.5; % center circle radius
Hlabel = Hbins;
nAngles = 16;
Dbins = linspace(0,360,nAngles+1);
nH = length(Hbins)-1;
nD = length(Dbins)-1;
Dbinc = mod(-90-Dbins,360) - 180;
id = find(diff(Dbinc)~=(Dbinc(2)-Dbinc(1)));
Dbinc(id+1:end) = Dbinc(id+1:end)-360; 

np = 20; % # of points for curve on polygons
X = nan*ones(np*2+1,nH*nD); % places holders for polygons
Y = X;
n = 1; % polygon counter
I = nan*ones(nD,nH); % empty intensity matrix
for jj=1:nH
   for kk=1:nD
      %if jj==1 % Build polygons.
         tmp = linspace(Dbinc(kk),Dbinc(kk+1),np)*pi/180;
         X(:,n) = [(Hbins(jj)+Hshift)*cos(tmp(1)) ...
            (Hbins(jj+1)+Hshift)*cos(tmp) ...
            (Hbins(jj)+Hshift)*cos(fliplr(tmp))];
         Y(:,n) = [(Hbins(jj)+Hshift)*sin(tmp(1)) ...
            (Hbins(jj+1)+Hshift)*sin(tmp)...
            (Hbins(jj)+Hshift)*sin(fliplr(tmp))];
         n = n + 1;
      %end
      if jj==nH
         id = find(B.H>=Hbins(jj) & B.D>=Dbins(kk) & ...
            B.D<Dbins(kk+1));
      else
         id = find(B.H>=Hbins(jj) & B.H<Hbins(jj+1) & ...
            B.D>=Dbins(kk) & B.D<Dbins(kk+1));
      end
      I(kk,jj) = length(id);
   end
end
I(I==0) = nan;

figure; orient landscape
hold on

fill(X,Y,I(:)'/N*100,'EdgeColor',[.8 .8 .8])
axis equal;
axis([-Hbins(end)-Hshift-.05 Hbins(end)+Hshift+.05...
   -Hbins(end)-Hshift-.05 Hbins(end)+Hshift+.05])
set(gca,'Visible','off')

[Hpx,Hpy] = pol2cart(45*pi/180,Hlabel+Hshift);
for jj=1:length(Hpx)
   text(Hpx(jj)+.03,Hpy(jj)+.07,[num2str(Hlabel(jj)),'m'],...
      'Color','k','FontSize',8,'FontWeight','bold',...
      'BackgroundColor','w');
end
h = plot(Hpx,Hpy,'k.','MarkerSize',8);
caxis([0 7.5])

hcb = colorbar('EastOutside');
p = get(hcb,'Position');
set(hcb,'Position',[p(1)+.05 p(2)+.1 .015 p(4)-.2],'FontSize',8)
set(get(hcb,'YLabel'),'String','% occurence','FontSize',10)

set(gca,'DefaultTextUnits','normalized')
title = ['Buoy Number: ' waves.buoy '; ' datestr(waves.stdmet.time(1), 'mm/dd/yyyy') ' to ' datestr(waves.stdmet.time(end), 'mm/dd/yyyy') ''];
text(.5,1.02,title,'FontSize',16,...
   'HorizontalAlignment','center')
%% Get the means in each bin.
%Dbinso = 90:Dbins(2):270;
Dbinso = 0:Dbins(2):360;
n = 1; 
N = length(B(1).time);
for ii=1:length(Hbins)-1
   for jj=1:length(Dbinso)-1
      if ii==length(Hbins)-1
         id = find(B.H>=Hbins(ii) & B.D>=Dbinso(jj) & ...
            B.D<Dbinso(jj+1));
      else
         id = find(B.H>=Hbins(ii) & B.H<Hbins(ii+1) & ...
            B.D>=Dbinso(jj) & B.D<Dbinso(jj+1));
      end
      statsHs(n).Hbin = Hbins(ii:ii+1);
      statsHs(n).Dbin = Dbinso(jj:jj+1);
      statsHs(n).percent = length(id)/N*100;
      statsHs(n).num = length(id);
      statsHs(n).H = nanmean(B.H(id));
      statsHs(n).Hstd = nanstd(B.H(id));
      statsHs(n).Hmax=max(B.H(id));
      statsHs(n).Hmin=min(B.H(id));
      statsHs(n).D = nanmean(B.D(id));
      statsHs(n).Dstd = nanstd(B.D(id));
      statsHs(n).T = nanmean(B.T(id));
      statsHs(n).Tstd = nanstd(B.T(id));
      statsHs(n).Tmax = max(B.T(id));
      statsHs(n).Tmin = min(B.T(id));
      statsHs(n).ws = nanmean(B.ws(id));
      statsHs(n).wsstd = nanstd(B.ws(id));
      statsHs(n).wd = nanmean(B.wd(id));
      statsHs(n).wdstd = nanstd(B.wd(id));
      n = n + 1;
      save("statsHs")
   end
end

outname = ['wave_climatology_' waves.buoy '.mat'];
save(outname, 'statsHs')
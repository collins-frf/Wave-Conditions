% plot results of the energy flux wave binning method

% call efm
EFM = efm(E);
centerRadius=0

%% Build bin lines
Dbins = naut2cart(EFM.Dmlims)-180;
Hbins = EFM.Hslims;
Hbins = [0*ones(1,size(Hbins,2)); Hbins];

I = EFM.binFreq(:);
I = I./sum(I)*100;

np = 20; % # of points for curve on polygons
nD = length(Dbins)-1;
nH = size(Hbins,1)-1;
X = nan*ones(np*2,nH*nD);
Y = X;
n = 1; % polygon counter
for dd=1:nD
   tmp = linspace(Dbins(dd),Dbins(dd+1),np)'*pi/180;
   for hh=1:nH
      a = [tmp;flipud(tmp)];
      r = [repmat(Hbins(hh,dd),[np,1]);repmat(Hbins(hh+1,dd),[np,1])];
      [X(:,n),Y(:,n)] = pol2cart(a,r);
      n = n + 1;
   end
end
[x,y] = pol2cart((naut2cart(E.Dm)-180)*pi/180,E.Hs+centerRadius);
%%
figureSize(10,7.5);
ax = axes;
hold on

% Plot the Bins
plot(x,y,'.','Color',[0.5 0.5 0.5]);
fill(X,Y,I,'EdgeColor','w','LineWidth',0.1,...
   'FaceAlpha',0.7);
axis equal;

% Set max Hs here!
Hbm = max(Hbins(:));
axis([-Hbm-centerRadius-.05 Hbm+centerRadius+.05...
   -Hbm-centerRadius-.05 0])
caxis([0 ceil(max(I))])

% Hs label
for jj=0.5:0.5:round(Hbm*2)/2
   xtmp = (jj+centerRadius).*cosd(naut2cart(90:270)-180);
   ytmp = (jj+centerRadius).*sind(naut2cart(90:270)-180);
   plot(xtmp,ytmp,'k','LineWidth',0.5,'LineStyle',':',...
      'Color',[0.5 0.5 0.5]);
   text(xtmp(1),ytmp(1),num2str(jj,'%.1f'),...
      'HorizontalAlignment','center','VerticalAlignment','bottom',...
      'FontSize',8)
end
text((Hbm+centerRadius+.05)/2,.15,'H_s (m)',...
   'HorizontalAlignment','center','VerticalAlignment','bottom',...
   'FontSize',8)

% Bin centers and Hs values
[Hpx,Hpy] = pol2cart((pi/180)*(naut2cart(EFM.binDm(:))-180),...
   EFM.binHs(:)+centerRadius);
plot(Hpx,Hpy,'ko','MarkerSize',5);
text(Hpx+.04,Hpy+.03,num2str(EFM.binHs(:),'%.2f'),...
   'FontWeight','bold','FontSize',8);

%h = text(0,0.35,[datestr(E.t(1),'mm/yyyy'),' - ',...
%   datestr(E.t(end),'mm/yyyy')],...
%   'FontSize',12,'HorizontalAlignment','center','FontWeight','bold');
h = text(0,0.35,name{:},...
   'FontSize',12,'HorizontalAlignment','center','FontWeight','bold');


% Colorbar
hcb = colorbar('SouthOutside');
p = hcb.Position;
hcb.Position = [p(1) 0.18 p(3) 0.02];
hcb.Label.String = '% occurence';
hcb.FontSize = 10;
cdata = hcb.Face.Texture.CData;
alphaVal = 0.7;
cdata(end,:) = uint8(alphaVal * cdata(end,:));
hcb.Face.Texture.ColorType = 'truecoloralpha';
hcb.Face.Texture.CData = cdata;

ax.Visible = 'off';

%%
print(['bins_',nametrim],'-dpng','-r300')

%%
save(['bins_',nametrim],'EFM');

%%
Hs = EFM.binHs';
Hs = Hs(:);
Tp = EFM.binTp';
Tp = Tp(:);
Dm = EFM.binDm';
Dm = Dm(:);
Occurrence = EFM.binFreq';
Occurrence = Occurrence(:)/sum(Occurrence(:));
bins = num2cell([repmat('bin',length(Hs),1),num2str((1:length(Hs))')],2);
dataTable = table(Hs,Tp,Dm,Occurrence,'RowNames',bins);

%writetable(dataTable,'efmBins','fileType','spreadsheet',...
%  'Sheet',[datestr(E.t(1),'mmyyyy'),'_',...
%  datestr(E.t(2),'mmyyyy')],'WriteRowNames',1)

writetable(dataTable,'efmBins','fileType','spreadsheet',...
  'Sheet',name{:},'WriteRowNames',1)
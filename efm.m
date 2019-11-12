function out = efm(E)
angE = 0;
angW = 360;

id = E.Dm>angE & E.Dm<angW;
Ec.angE = angE;
Ec.angW = angW;
Ec.t = E.t(id);
Ec.Hs = E.Hs(id);
Ec.Tm = E.Tm(id);
Ec.Tp = E.Tp(id);
Ec.Dm = E.Dm(id);
Ec.Co = (9.81/(2*pi))*Ec.Tp;
Ec.Cg = 0.5*Ec.Co;
Ec.Ef = (1025*9.81*Ec.Hs.^2/8).*Ec.Cg;
%%
nD = 3; % # of direction bins
nH = 3; % # of height bins
%nD = 5; % # of direction bins
%nH = 4; % # of height bins

%% Determine equal energy direction bins.
binED = nan*ones(nD,1); % % of energy in direction bin... should all 
                        % equal 0.2 for na=5
Dmlims = nan*ones(nD,1);  % angles that define the bins
for ii=1:nD
   % Find the points between current east angle and angW
   if ii==1
      D1 = Ec.angE;
   else
      D1 = Dmlims(ii-1);
   end
   da = Ec.angW - D1;
   D2 = D1 + da;
   id = Ec.Dm>=D1 & Ec.Dm<D2;
   bE = sum(Ec.Ef(id))/sum(Ec.Ef); % % of energy in this current bin
   
   % Continue to narrow the bin until the energy in it is 1/nD
   if ii<nD
      while bE>(1/nD)
         da = da - 0.1;
         D2 = D1 + da;
         id = Ec.Dm>=D1 & Ec.Dm<D2;
         bE = sum(Ec.Ef(id))/sum(Ec.Ef);
         %if bE<(1/nD)+0.015
      end
   end
   
   binED(ii) = bE;   % Again, these should all be 0.2 for nD=5
   Dmlims(ii) = D2;  % angles that define the bins
   %plot([ang2 ang2],[0 7],'k')
end
Dmlims = [Ec.angE; Dmlims];

%% Now, within each direction bin, make equal energy wave height bins.
binEH = nan*ones(nH,nD);   % % of energy in height bin... should all be 
                           % equal to 1/(na*nH)
Hslims = nan*ones(nH,nD);  % heights that define the bins
binHs = Hslims;
binTp = Hslims;
binDm = Hslims;
binFreq = Hslims;
binlims = nan*ones(5,2,nD,nH);
nN = 12;%nD*nH;
for ii=1:nD
   for jj=1:nH
      
      % Find the points between current H and max H in direction bin
      if jj==1
         H1 = 0;
         %H1 = 0.5;
      else
         H1 = Hslims(jj-1,ii);
      end      
      D1 = Dmlims(ii);
      D2 = Dmlims(ii+1);
      id = Ec.Dm>D1 & Ec.Dm<D2;
      dH = max(Ec.Hs(id)) - H1;
      H2 = H1 + dH;
      
      if jj==nH
         id = Ec.Dm>D1 & Ec.Dm<D2 & Ec.Hs>=H1 & Ec.Hs<=H2;
      else
         id = Ec.Dm>D1 & Ec.Dm<D2 & Ec.Hs>=H1 & Ec.Hs<H2;
      end
      bE = sum(Ec.Ef(id))/sum(Ec.Ef);
      
      if bE==0
         %nN = nN - 1;
         continue
      end
      % Continue to narrow the bin until the energy in it is 1/(nD*nH)
      if jj<nH
         while bE>1/(nD*nH)
            dH = dH - 0.01;
            H2 = H1 + dH;
            id = Ec.Dm>=D1 & Ec.Dm<D2 & Ec.Hs>=H1 & Ec.Hs<H2;
            bE = sum(Ec.Ef(id))/sum(Ec.Ef);
         end
      end
      Ef = Ec.Ef(id);
      %Hs = Ec.Hs(id);
      Tp = Ec.Tp(id);
      Dm = Ec.Dm(id);
      
      % # of points in bin
      binFreq(jj,ii) = length(find(id));
      % mean Hs in bin
      binHs(jj,ii) = sqrt(mean(Ef)*32*pi/(1025*9.81^2*mean(Tp)));
      % mean Tp in bin
      binTp(jj,ii) = mean(Tp);
      % mean Dm in bin
      binDm(jj,ii) = mean(Dm);
      % Energy in bin
      binEH(jj,ii) = bE;
      % Hs that defins bin
      Hslims(jj,ii) = H2;
      
      binlims(:,:,ii,jj) = [[D1 D2 D2 D1 D1]',[H1 H1 H2 H2 H1]'];
      %plot([angs(ii) angs(ii+1) angs(ii+1) angs(ii) angs(ii)],...
      %   [H1 H1 H2 H2 H1],'Color',[0.5 0.5 0.5])
      %plot(binDm(jj,ii),binHs(jj,ii),'yo')
      %plot([angs(ii) angs(ii+1) angs(ii+1) angs(ii) angs(ii)],...
      %   [H1 H1 H2 H2 H1],'k')
      %plot(binDm(jj,ii),binHs(jj,ii),'ro')
   end
end

out.Hslims = Hslims;
out.Dmlims = Dmlims;
out.binlims = binlims;
out.binHs = binHs;
out.binDm = binDm;
out.binTp = binTp;
out.binFreq = binFreq;
out.binED = binED;
out.binEH = binEH;










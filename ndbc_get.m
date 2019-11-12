function out = ndbc_get(buoy,dates,dtype,method,dirres)
%
%NDBC_GET  Gets NOAA NDBC data from DODS server.
%
%   OUT = NDBC_GET(BUOY,DATES) gets the Standard Metorological data for
%   BUOY from (or between and including) DATE(S) where BUOY is an NDBC id
%   number and DATE(S) is/are Matlab datenum(s), a year, years or '*'
%   (get everything available).
%
%   OUT = ndbc_get(BUOY,DATES,DTYPE) gets the Standard Meteorological
%   data for DTYPE='smet' (default), gets the Spectral data for
%   DTYPE='spec' or gets both for DTYPE='all'.
%
%   OUT = ndbc_get(...,METHOD,DIRRES) for METHOD='lh' (default) uses the
%   Longuet-Higgins et al. (1963) method for calculating the directional
%   distribution function. METHOD='lhm' is a weighted version of the LH
%   method to avoid negative energies. METHOD='MEM' uses the Maximum
%   Entropy Method of Lygre & Krogstad (1986 - JPO). METHOD='MLM' uses the
%   Maximum Likelihood Method of Oltman-Shay & Guza (1984 - JPO). DIRRES
%   (default=1) is the direction resolution to use for any given method.
%
%   Input:
%     BUOY = num (ex. buoy=41010)
%     DATE(S) = datenum | [datenum1 datenum2] | year | [year1 yearend] | '*';
%     DTYPE = 'smet' (default) | 'spec' | 'all'
%     METHOD = 'lh' (default) | 'lhm' | 'mem' | 'mlm'
%     DIRRES = num (default=1)
%   Output:
%     OUT = structure:
%           buoy = BUOY
%           longitude = longitude
%           latitude = latitude
%           depth = depth in m
%           comment = comment in netcdf file (at bottom)
%           stdmet = structure: Standard Meteorological Data
%              See http://www.ndbc.noaa.gov/measdes.shtml for descriptions.
%           swden = structure: Spectral wave density data
%
%   NOTE: MWD is in Nautical convention: direction from measured cw from
%   North.
%
%   Requires:
%     nctoolbox
%
% Dave Thompson (dthompson@usgs.gov)
%

%% Check inputs
more off
buoy = lower(num2str(buoy));

if ~exist('dtype','var')
   dtype = 'smet';
end
data = {'stdmet','swden';'h','w'};
if strcmp('smet',dtype)
   data = data(:,1);
elseif strcmp('spec',dtype)
   data = data(:,2);
end

if ~exist('method','var')
   method = 'lh';
end
if ~exist('dirres','var')
   dirres = 1;
end
dirs = 0:dirres:359;

%% Year.
out.buoy = buoy;
out.comment = [];
out.longitude = [];
out.latitude = [];
out.depth = [];

%if ~exist('dates','var')
tmp = webread(['http://dods.ndbc.noaa.gov/thredds/catalog/data/stdmet/',...
   buoy,'/catalog.html']);
   disp(' Reading location and comment.')
   tmp = regexp(tmp,[buoy,'\w*\.\w+'],'match');
   tmp = unique(tmp);
   nc = ncgeodataset(['https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/',...
      buoy,'/',tmp{1}]);
   out.comment = nc.attribute('comment');
   out.longitude = nc.data('longitude');
   out.latitude = nc.data('latitude');
%end

% Get the water depth
disp(' Reading depth.')
url = ['http://www.ndbc.noaa.gov/station_page.php?station=',...
   num2str(buoy)];
tmp = urlread(url);
dep = regexp(tmp,'Water depth:[</b>]* \d+\.*\d*','match');
if ~isempty(dep)
   dep = regexp(dep{:},'\d+\.*\d*','match');
   out.depth = str2num(dep{:});
else
   out.depth = nan;
end

if ~exist('dates','var')
   return
end

dvec = datevec(now);
y = dvec(1);

if dates=='*' % Get all data.
   url = ['http://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/',...
      buoy,'/catalog.xml'];
   tmp = urlread(url);
   tmp = regexp(tmp,[buoy,'h\d{4}\.nc'],'match');
   tmp = unique(tmp);
   tmp = regexp(tmp,'h\d{4}','match');
   tmp = regexp([tmp{:}],'\d{4}','match');
   tmp = [tmp{:}];
   years = str2num(vertcat(tmp{:}));
   if years(end)==9999
      tmp = years(end-1);
   else
      tmp = years(end);
   end
   disp(' Data for:');
   fprintf('\t%.0f\n',years);
   dates = [datenum(years(1),1,1) datenum(tmp,12,31)];
   
else % Use specific dates that came in.
   if length(dates)==1
      if length(num2str(dates))==4 % Get one year.
         dates = [datenum(dates,1,1) datenum(dates,12,31,23,59,0)];
      else % Get one day.
         dates = [dates dates+1];
      end
   else
      dates = dates([1 end]);
   end
   
   if length(num2str(dates(1)))==4
      years = (dates(1):dates(2))';
      dates = [datenum(dates(1),1,1) datenum(dates(2)+1,1,1)];
   else
      dvecin = datevec(dates);
      years = unique(dvecin(:,1));
      years = (years(1):years(end))';
      
      if dates(1)>(now-45)
         years = 9999;
      elseif dates(2)<(now-45)
         %years = years;
      else
         tmp = datevec(datestr(now-45));
         if (tmp(2)==11 || tmp(2)==12) && diff(dvecin(:,1))==1
            years = years(1:end-1);
         end
         years = [years;9999];
      end
   end
end

%       tmp = datevec(dates);
%       years = unique(tmp(:,1));
%       years = (years(1):years(end))';
%    end
%    url = ['http://dods.ndbc.noaa.gov/thredds/dodsC/data/',...
%       data{1,1},'/',buoy,'/',buoy,data{2,1},num2str(years(end)),'.nc'];
%    nc = ncgeodataset(url);
%    if ~isempty(nc)
%       t = epoch2datenum(double(nc.data('time')));
%       if t(end)<dates(1)
%          years = 9999;
%       elseif t(end)>dates(2)
%          years = years;
%       else
%          years = [years;9999];
%       end
%    end
% end
% end

%% Get the data.

for jj=1:size(data,2)
   out.(data{1,jj}).time = [];
   if strcmp('swden',data{1,jj})
      out.(data{1,jj}).frequency = [];
      out.(data{1,jj}).direction = [];
   end
end

for ii=1:length(years)
   fprintf(['\tGetting data from ',num2str(years(ii)),'...\n'])
   for jj=1:size(data,2)
      if years(ii)~=9999
         fprintf(['\t\t',data{1,jj},'\n'])
         url = ['http://dods.ndbc.noaa.gov/thredds/dodsC/data/',...
            data{1,jj},'/',buoy,'/',buoy,data{2,jj},num2str(years(ii)),'.nc'];
         
         % Open url.
         nc = ncgeodataset(url);
         
         % Get time.
         t = epoch2datenum(double(nc.data('time')));
         
         % Get variables.
         vars = nc.variables;
         if strcmp('stdmet',data{1,jj})
            vars = setxor(vars,{'time','longitude','latitude'});
         else
            vars = setxor(vars,{'time','longitude','latitude','frequency'});
         end
         
         % Need a catch if it is currently January and you are trying to
         % get December of last year. I've found that the time vector
         % stored in the netCDF file will include the entire year, but the
         % data does not.
         dend = datevec(dates(end));
         dnow = datevec(now);
         if dend(1)==dnow(1) && dend(2)==1 && dend(3)==1 && dnow(2)==1
            tmp = nc.data(vars{1});
            if length(tmp)~=length(t)
               t = t(1:length(tmp));
            end
         end
         
         if ii==1
            if length(dates)==1
               idt = length(t);
            else
               idt = find(t>=dates(1) & t<=dates(2));
            end
         else
            idt = find(t>=(max([dates(1); out.stdmet.time]))...
               & t<=dates(2));
         end
         
         % Initialize directional stuff
         if ii==1 && strcmp('swden',data{1,jj})
            out.(data{1,jj}).frequency = nc.data('frequency');
            out.swden.direction = dirs';
         end
         
         % Note that there is a frequency discrepancy between historical data
         % and recent (9999) data... historical has 47 fs while recent has 46.
         % On top of that the first 13 fs of the recent data are different
         % than the historical in the thousanths place. My solution... pad
         % recent data with nans and use historical fs.
         % Skip what is right here and take care of this down where the data
         % is read in.
         
         

         if ~isempty(idt)
            % Store time.
            out.(data{1,jj}).time = [out.(data{1,jj}).time; t(idt)];
            
            % Get data from each variable.
            for kk=1:length(vars)
               % If the variable doesn't yet exist then initialize it.
               if ~isfield(out.(data{1,jj}),vars{kk})
                  out.(data{1,jj}).(vars{kk}) = [];
               end
               
               % axes for stdmet vars are time, lat, lon.
               % axes for swden vars are time, freq, lat, lon.
               % lat and lon are singular. This statement will work for
               % both stdmet (grabbing time, lat) and for swden (grabbing
               % time, freq)
               tmp = nc{vars{kk}}(idt,:);
               %if size(tmp,2)~=length(out.swden.frequency)
               %   tmp = [ones(size(tmp,1),1)*nan tmp];
               %end
               %end
               % Store it.
               out.(data{1,jj}).(vars{kk}) = ...
                  [out.(data{1,jj}).(vars{kk}); double(tmp)];
            end
         end
      else
         url = ['http://www.ndbc.noaa.gov/data/realtime2/',num2str(buoy),...
            '.txt'];
         cols = [6:17 19];   % Data columns.
         tmp = webread(url);
         head = tmp(1:regexp(tmp,'\d','once')-1);
         head = regexp(head,'\n','split');
         head = regexp(head{1},'#*\w+','match');
         vars = head(cols);
         
         dat = tmp(regexp(tmp,'\d','once'):end);
         dat = strrep(dat,'MM','nan'); % This is for 'recent data'.
         dat = textscan(dat,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f\n','Delimiter',' ','MultipleDelimsAsOne',1);
         t = datenum([dat{1} dat{2} dat{3} dat{4} dat{5} zeros(length(dat{1}),1)]);
         dat = dat(cols);
         
         vars = [vars;{'wind_dir','wind_spd','gust',...
               'wave_height','dominant_wpd','average_wpd','mean_wave_dir',...
               'air_pressure','air_temperature','sea_surface_temperature',...
               'dewpt_temperature','visibility','water_level'}];
         [~,id] = sort(vars(2,:));
         vars = vars(:,id);
         dat = dat(id);
         
         if ii==1
            idt = find(t>=dates(1) & t<=dates(2));
         else
            idt = find(t>=(max([dates(1); out.stdmet.time]))...
               & t<=dates(2));
         end
         
         if ~isempty(idt)
            % Store time.
            out.(data{1,jj}).time = [out.(data{1,jj}).time; flipud(t(idt))];
            % Get data from each variable.
            for kk=1:length(vars)
               % If the variable doesn't yet exist then initialize it.
               if ~isfield(out.(data{1,jj}),vars{2,kk})
                  out.(data{1,jj}).(vars{2,kk}) = [];
               end
               out.(data{1,jj}).(vars{2,kk}) = ...
                  [out.(data{1,jj}).(vars{2,kk}); flipud(dat{kk}(idt))];
            end
         end
      end
   end
end

%%
if exist('progressBarAscii')==2
   pba = 1;
else
    pba = 0;
end

if sum(sum(strcmp('swden',data)))
   % Calculate the directional distribution function.
   out.swden.ddf = ones([size(out.swden.mean_wave_dir)...
      length(dirs)])*nan;
   fprintf('\tCalculating directional distribution function...\n')
   switch method
      case {'lh'}
         if pba==1; pbar = progressBarAscii(length(dirs)); end
         for mm=1:length(dirs)
            out.swden.ddf(:,:,mm) = dirres*(1/180)*(1/2 + ...
               out.swden.wave_spectrum_r1.*cosd(dirs(mm) - ...
               out.swden.mean_wave_dir) + ...
               out.swden.wave_spectrum_r2.*cosd(2*(dirs(mm) - ...
               out.swden.principal_wave_dir)));
            if exist('pbar','var'); pbar(mm); end
         end
      case {'lhm'} % LHM with weighting to aviod negative energies.
         if pba==1; pbar = progressBarAscii(length(dirs)); end
         for mm=1:length(dirs)
            out.swden.ddf(:,:,mm) = dirres*(1/180)*(1/2 + ...
               (2/3)*(out.swden.wave_spectrum_r1).*cosd(dirs(mm) - ...
               out.swden.mean_wave_dir) + ...
               (1/6)*(out.swden.wave_spectrum_r2).*cosd(2*(dirs(mm) - ...
               out.swden.principal_wave_dir)));
            if exist('pbar','var'); pbar(mm); end
         end
      case {'mlm'} % MLM - This routine needs to be vectorized.
         warning(' mlm takes a LOOOOOONG time!')
         if pba==1
            pbar = progressBarAscii(size(out.swden.ddf,1)*...
               size(out.swden.ddf,2));
         end
         zz = 1;
         for kk=1:size(out.swden.ddf,1)
            for ll=1:size(out.swden.ddf,2)
               out.swden.ddf(kk,ll,:) = (1/dirres)*...
                  spec_mlm(out.swden.wave_spectrum_r1(kk,ll),...
                  out.swden.wave_spectrum_r2(kk,ll),...
                  out.swden.mean_wave_dir(kk,ll),...
                  out.swden.principal_wave_dir(kk,ll),dirs);
               if exist('pbar','var'); pbar(zz); end
               zz = zz + 1;
            end
         end
      case {'mem'}  % MEM
         tmp = spec_mem(out.swden.wave_spectrum_r1/100,...
            out.swden.wave_spectrum_r2/100,...
            out.swden.mean_wave_dir,out.swden.principal_wave_dir,dirs);
         out.swden.ddf = (1/dirres)*permute(tmp,[3 1 2]);
   end
   out.swden.method = method;
   
   % Finally calculate the 2D spectra.
   %out.swden.t = out.swden.time;
   out.swden.S = permute(repmat(out.swden.spectral_wave_density,...
      [1 1 size(out.swden.ddf,3)]).*out.swden.ddf,...
      [2 3 1]);
   
   % Calculate Hs.
   % Have to force Nans to zero to use trapz.
   tmp = out.swden.S;
   tmp(isnan(tmp)) = 0;
   out.swden.Hs = squeeze(4*sqrt(trapz(trapz(out.swden.frequency,...
      tmp)*dirres)));
   out.swden.Hs(squeeze(sum(sum(tmp)))==0) = nan;
   
   tmp = reshape(out.swden.S,...
      [size(out.swden.S,1)*size(out.swden.S,2) size(out.swden.S,3)]);
   [tmp,id] = max(tmp);
   [i,j] = ind2sub([size(out.swden.S,1) size(out.swden.S,2)],id);
   out.swden.Tp = (1/out.swden.frequency(i))';
   out.swden.Tp(isnan(tmp)) = nan;
   out.swden.Dp = out.swden.direction(j);
   out.swden.Dp(isnan(tmp)) = nan;
   
end



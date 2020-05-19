function [wave] = getwaveFRF(d1, d2, gnum)
% %     function takes gauge number and returns data from FDIF server quick version in matlab based on that done in python
%   This function grabs data from the THREDDS server at CHL(2) for waves.
%   This code is meant to be used as an example, throurgh debugging has not been done
%   for this set of scripts.  The more complete version exists in python for the coastal model 
%   Test Bed (CMTB) 
%   Written by: Spicer Bak, PhD
%   email: Spicer.bak@usace.army.mil
%   edits by Julia fiedler, UCSD
%   
%   
%     INPUTS
%       d1=datenum(2015,10,3); example
%       d2=datenum(2015,10,4); example
%       svrloc: 
%       gnum:  gauge number of interset 
%         1 = waverider 430  - 26 m
%         2 = waverider 630  - 17 m 
%         if gnum==1 %26 m wavericder
%            gnum==2 % 17 m waverider
%            gnum==3 % pressure sensor at x=200m
%            gnum==4 % pressure sensor at x=150m
%            gnum==5 % pressure sensor at x=125m
%            gnum==6 % pressure sensor at x=100m
%             gnum==7 % 6 m awac
%             gnum==8 % 4.5 m awac
%             gnum==9 % 11 m awac
%             gnum==10 % 8 m array
%         http://chldata.erdc.dren.mil/thredds/catalog/frf/catalog.html 
%         to find other gauges of interest
% 
%     RETURNS
%        data structure data with data inside

%% Matlab Get Data function start
if d2<d1
    disp(' your times are backwards')
    return 
end
%% set url's 
%svrloc='http://134.164.129.55/thredds/dodsC/FRF/';  % The prefix for the CHL thredds server
svrloc='https://chldata.erdc.dren.mil/thredds/dodsC/frf/';  % The prefix for the CHL thredds server

% add other wave gauges here
% TODO: change these arbitrary numbers to actual string calls
if gnum==1
    urlback='oceanography/waves/waverider-26m/waverider-26m.ncml'; %26 m wavericder
elseif gnum==2
    urlback='oceanography/waves/waverider630/waverider-17m.ncml'; % 17 m waverider
elseif gnum==3
    urlback='oceanography/waves/xp200m/xp200m.ncml'; % pressure sensor at x=200m
elseif gnum==4
    urlback='oceanography/waves/xp150m/xp150m.ncml'; % pressure sensor at x=150m
elseif gnum==5
    urlback='oceanography/waves/xp125m/xp125m.ncml'; % pressure sensor at x=125m
elseif gnum==6
    urlback='oceanography/waves/xp100m/xp100m.ncml'; % pressure sensor at x=100m
elseif gnum==7
    urlback='oceanography/waves/awac-6m/awac-6m.ncml'; % 6 m awac
elseif gnum==8
    urlback='oceanography/waves/awac-4.5m/awac-4.5m.ncml'; % 4.5 m awac
elseif gnum==9
    urlback='oceanography/waves/awac-11m/awac-11m.ncml'; % 11 m awac
elseif gnum==10
    urlback='oceanography/waves/8m-array/8m-array.ncml'; % 8 m array
else
    disp ' go to http://chldata.erdc.dren.mil/thredds/catalog/frf/catalog.html and browse to the gauge of interest and select the openDAP link'
end
url=strcat(svrloc,urlback);
%% Main program
time=ncread(url,'time'); % downloading time from server
%tunit= ncreadatt(url,'time','units'); % reading attributes of variable time  # how to read an attribute
% converting time to matlab datetime
mtime=time/(3600.0*24)+datenum(1970,1,1);
% fprintf('Parsing time\nWave Record at this location starts: %s  ends: %s\n',datestr(min(mtime)),datestr(max(mtime)))
% finding index that corresponds to dates of interest
itime=find(d1 <= mtime & d2>= mtime); % indicies in netCDF record of data of interest
% disp('Retrieving Data')
% pulling data that is time length dependent (D1,D2)
% TODO: fix these try/catch statements
if ~isempty(itime)
    try
        wave.time=mtime(itime);  % record of wave time indicies in matlab datetime format
        wave.Hs=ncread(url,'waveHs',min(itime),length(itime));
        wave.depth=ncread(url,'depth',min(itime),length(itime));
        wave.station_name=ncreadatt(url,'/','title');
    catch
    end
    
    try
        wave.lat=ncread(url,'lat');
        wave.lon=ncread(url,'lon');
    catch
        wave.lat=ncread(url,'latitude');
        wave.lon=ncread(url,'longitude'); 
    end
    
    try
        wave.HsIG = ncread(url,'waveHsIg',min(itime),length(itime));
    catch
    end
    
    try
        wave.fp=ncread(url,'waveFp',min(itime),length(itime));
    catch
        wave.fp=1/ncread(url,'waveTp',min(itime),length(itime));
    end
    
    try
        wave.Tp=ncread(url,'waveTp',min(itime),length(itime));
    catch
        wave.Tp=1/ncread(url,'waveTp',min(itime),length(itime));
    end
    
    try
        wave.Dp=ncread(url,'wavePeakDirectionPeakFrequency',min(itime),length(itime));
    catch
    end
    
    try
        wave.spec2D =permute(ncread(url, 'directionalWaveEnergyDensity' , [1,1,min(itime)],[inf,inf,length(itime)]),[3, 1, 2]); % arranging with hours in first index
        if length(itime) ==1
            wave.spec2D = squeeze(wave.spec2D);
        end
    catch
    end
    
    try
        wave.dirpeak=ncread(url,'waveDp',min(itime),length(itime));
    catch
    end
    
    try
        wave.dirpeak=ncread(url,'wavePeakDirectionPeakFrequency',min(itime),length(itime));
    catch
    end
    
    try
        wave.spec1D=permute(ncread(url,'waveEnergyDensity',[1,min(itime)],[inf,length(itime)]),[2, 1]); % arraging with hours in first value
    catch
    end
 
    try
        wave.freqSpread=ncread(url,'spectralWidthParameter',min(itime),length(itime)); % arraging with hours in first value
    catch
    end
    
    try
        wave.dirSpread=ncread(url,'directionalPeakSpread',min(itime),length(itime)); % arraging with hours in first value
    catch
    end
    
    try
        %non time scale dependent variables
        wave.frqbin=ncread(url,'waveFrequency');
        wave.dirbin=ncread(url,'waveDirectionBins');
    catch
    end
    
    %disp 'Data successfully grabbed'
    
else
    wave.time=0;
    wave.error = 'no wave data';
    fprintf('ERROR: There''s no wave data at %s during %s to %s\nTry another gauge\n', ncreadatt(url,'/','title'), datestr(d1), datestr(d2))
end



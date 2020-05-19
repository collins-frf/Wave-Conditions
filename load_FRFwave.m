function load_FRFwave(fname)
%% loads a flat netcdf file or a data structure from get data frf 
try  % loads from flat file
    waveTime=ncread(fname,'time');
    Hs=ncread(fname,'waveHs');
    Tp=ncread(fname,'waveTp');
    Dp=ncread(fname,'waveMeanDirection');
    nominaldepth=ncread(fname,'nominalDepth');
    wavedepth=ncread(fname,'depth');
    waveFrequency=ncread(fname,'waveFrequency');
    waveDirection=ncread(fname,'waveDirectionBins');
    waveEnergyDensity=ncread(fname,'directionalWaveEnergyDensity');
catch  % this will parse a structure from the getWaveFRF, from thredds
    waveTime = fname.time;
    Hs = fname.Hs;
    Tp = fname.Tp;
    Dp = fname.Dp;
    nominaldepth = median(fname.depth);  % incase there are multiple values
    wavedepth = fname.depth;
    waveFrequency = fname.frqbin;
    waveDirection = fname.dirbin;          
    waveEnergyDensity = fname.spec1D;  % this likely has to be permuted back to time last then netCDF convention for freq and dir
end
eval(['save FRFwave_forecast.mat waveTime Hs Tp Dp nominaldepth wavedepth waveFrequency waveDirection waveEnergyDensity'])

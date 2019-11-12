function dnum = epoch2datenum(epoch)
% epoch2datenum.m  A function to convert Epoch (Unix) time stamps to Matlab
%                  date numbers.
%
% USAGE:        dnum = epoch2datenum(epoch);
%
%       WHERE:  epoch is a 1D array of Epoch (Unix) time stamps (seconds
%                     since January 1, 1970 at midnight).
%               dnum is a 1D array of Matlab date numbers (days since
%                    January 1, 0000 at midnight).
%
% HISTORY:
% C. Sullivan,  09/10/07,   version 1.0

epochTimeBase = datenum(1970,1,1,0,0,0);
dnum = (epoch/(3600*24)) + epochTimeBase;
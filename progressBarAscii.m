function f = progressBarAscii(nMax)
%progressBar    Ascii progress bar.
%   progBar = progressBar(nSteps) creates a progress bar and returns a 
%   pointer to a function handle which can then be called to update it.
%
%   To update, call progBar(currentStep)
%
%   Example:
%      n = 5000; 
%      progBar = progressBar(n);
%      for tmp = 1:n
%        progBar(tmp); 
%      end
 
%   by David Szotten 2008
%   $Revision: 4 $  $Date: 2011-01-20 12:51:13 -0500 (Thu, 20 Jan 2011) $

lastPercentileWritten = 0;

fprintf('| 0%%');
for tmp = 1:90
	fprintf(' ');
end
fprintf('100%% |\n');

f = @updateBar;
	function updateBar(nCurrent)
		
		%what percentile are we up to
		currentPercentile = round(100*nCurrent/nMax);

		%have we passed another percentile?
		if (currentPercentile > lastPercentileWritten )
			
			%we may have increased by several percentiles,
			%so keep writing until we catch up
			percentileToWrite = lastPercentileWritten + 1;
			while(percentileToWrite <= currentPercentile)

				%for every 10th, use a '+' instead
				if mod(percentileToWrite,10)==0 || percentileToWrite==1
					fprintf('+');
					
					%are we done?
					if percentileToWrite == 100
						fprintf('\n');
					end
				else
					fprintf('.');
				end
				percentileToWrite = percentileToWrite + 1;
			end
			
			%update status
			lastPercentileWritten = currentPercentile;
		end
	end

end

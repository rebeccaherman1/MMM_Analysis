function [periods, mean_periodogram, low, high, ml, mh, sds] = ps_multi(hi, ti, c, PAD, groupings, groupings2)
%POWERSPECTRUM hi = sequence, a COLUMN vector, or the COLUMNS of a matrix. ti = time
%returns periodogram and frequencies. also plots periodogram vs.
%frequencies. Use add_to_fig=false to add to previously existing figure.
%should come in with a mean of 0.
grouping = true;
if(nargin<=5)
    groupings = nan; groupings2=nan;
    grouping = false;
end
sz = size(hi); N_scale = sum(~isnan(hi),1);
hi(isnan(hi))=0;
%hi_pad = hi_pad(1:min(sum(~isnan(hi_pad),1)), :);%hi(~isnan(hi(:,1)),:); %TODO don't throw away data...
hi_pad = [hi;zeros(PAD-sz(1),sz(2))];
%length(hi_pad(:,1));
N = PAD; 
ft = fft(hi_pad); %uses columns!!!!
df = 1/((N-1)*(ti(2)-ti(1)));%(ti(end)-ti(1));
%fmax = df*N/2;
num_freqs = floor(N/2); % = fmax/df;
freqs = [0:num_freqs, -num_freqs+1:-1]'*df;
positive_frequencies = freqs(2:num_freqs); periods = 1./positive_frequencies;
periodogram = ft.*conj(ft)./N_scale;%^2;%necessary to satisfy Parseval's theorem
positive_periodograms = periodogram(2:num_freqs,:);
if(grouping)
    model_mean_periodograms = splitapply(@(X) mean(X,2), positive_periodograms, groupings);
    institution_mean_periodograms = splitapply(@(X) mean(X,2), model_mean_periodograms, groupings2);
    mean_periodogram = mean(institution_mean_periodograms, 2);
else
    mean_periodogram = mean(positive_periodograms, 2);
    model_mean_periodograms = positive_periodograms;
end
[low, high] = confidence_interval(model_mean_periodograms');
[ml, mh] = confidence_interval(model_mean_periodograms',c);
sds = std(positive_periodograms,0,2)/sqrt(length(hi_pad(1,:)));
end


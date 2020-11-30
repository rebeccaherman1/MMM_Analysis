function [periods, mean_periodogram, low, high, ml, mh, sds] = ps_multi(seq, T, confidence, PAD, dim, groupings, groupings2)
%currently doesn't do any plotting. Is it OK to stick with this mean
%method for the simulated versions? 
grouping = true;
if(nargin<5)
    dim = 1;
end
if(nargin<6)
    groupings = nan; groupings2=nan;
    grouping = false;
end
seq = seq - mean(seq, dim);
N_scale = sum(~isnan(seq),dim);
seq(isnan(seq))=0; 

ft = fft(seq, PAD, dim);
df = 1/((PAD-1)*(T(2)-T(1)));
num_freqs = floor(PAD/2);
freqs = [0:num_freqs, -num_freqs+1:-1]*df;
positive_frequencies = freqs(2:num_freqs); periods = 1./positive_frequencies;
periodogram = ft.*conj(ft)./N_scale;%^2;%necessary to satisfy Parseval's theorem
switch dim
    case 1
        m_dim = 2;
        positive_periodograms = periodogram(2:num_freqs,:);
    otherwise
        m_dim = 1;
        positive_periodograms = periodogram(:,2:num_freqs);
end
%this is not a weighted mean because the PS don't lose meaningful variance
%as they are averaged. Of course, those with more simulations still have
%less noise...
if(grouping)
    model_mean_periodograms = splitapply(@(X) mean(X,m_dim), positive_periodograms, groupings);
    institution_mean_periodograms = splitapply(@(X) mean(X,m_dim), model_mean_periodograms, groupings2);
    mean_periodogram = mean(institution_mean_periodograms, m_dim);
else
    mean_periodogram = mean(positive_periodograms, m_dim);
    model_mean_periodograms = positive_periodograms;
end
[low, high] = confidence_interval(model_mean_periodograms, m_dim);
[ml, mh] = confidence_interval(model_mean_periodograms, m_dim, confidence);
sds = std(positive_periodograms,0,2)/sqrt(size(seq, dim));
end


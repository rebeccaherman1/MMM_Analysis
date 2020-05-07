function [ft, freqs, positive_periodogram, positive_frequencies ] = ps_mean_sahel(hi, ti, name, color, PAD, lw)
%POWERSPECTRUM hi = sequence, a COLUMN vector, or the COLUMNS of a matrix. ti = time
%returns periodogram and frequencies. also plots periodogram vs.
%frequencies. Use add_to_fig=false to add to previously existing figure.
if(nargin < 7)
    lw = 1;
end
%MAX_PERIOD = 103; 
hi_anomaly = hi - nanmean(hi,1);
sz = size(hi_anomaly); N_scale = sum(~isnan(hi_anomaly),1);
hi_anomaly(isnan(hi_anomaly))=0;
hi_pad = [hi_anomaly;zeros(PAD-sz(1),sz(2))];
N = PAD;
ft = fft(hi_pad); %uses columns!!!!
df = 1/((N-1)*(ti(2)-ti(1)));%(ti(end)-ti(1));
%fmax = df*N/2;
num_freqs = floor(N/2); % = fmax/df;
freqs = [0:num_freqs, -num_freqs+1:-1]'*df;
positive_frequencies = freqs(2:num_freqs); periods = 1./positive_frequencies;
periodogram = ft.*conj(ft)./N_scale;%^2;%necessary to satisfy Parseval's theorem
positive_periodogram = periodogram(2:num_freqs,:);
positive_ft = ft(2:num_freqs,:);
mean_periodogram = mean(positive_periodogram, 2);
plot(periods, mean_periodogram, 'color', color, 'Linewidth', lw, 'DisplayName', name); 
hold on;
xlabel('Period (years)')
ylabel('Power')
end


function [ft, freqs, positive_periodogram, positive_frequencies ] = ps_sahel(hi, ti, name, style, color, lw)
%POWERSPECTRUM hi = sequence, a COLUMN vector, or the COLUMNS of a matrix. ti = time
%returns periodogram and frequencies. also plots periodogram vs.
%frequencies. Use add_to_fig=false to add to previously existing figure.
if(nargin < 5)
    lw = 1;
end
N = length(hi(:,1)); N_norm = length(ti(:,1));    
hi = hi(~isnan(hi(:,1)),:); 
ft = fft(hi); %uses columns!!!!
df = 1/((N-1)*(ti(2)-ti(1)));%(ti(end)-ti(1));
%fmax = df*N/2;
num_freqs = floor(N/2); % = fmax/df;
freqs = [0:num_freqs, -num_freqs+1:-1]'*df;
positive_frequencies = freqs(2:num_freqs); periods = 1./positive_frequencies;
periodogram = ft.*conj(ft)/N_norm;%^2;%necessary to satisfy Parseval's theorem
positive_periodogram = periodogram(2:num_freqs,:);
positive_ft = ft(2:num_freqs,:);
s = size(positive_periodogram);
for i = 1:s(2)
    if i==1
        plot(periods, positive_periodogram(:,i), style, 'color', color, 'Linewidth', lw, 'DisplayName', name); 
    else
        plot(periods, positive_periodogram(:,i), style, 'color', color, 'Linewidth', lw, 'HandleVisibility','off');
    end
    hold on;
end
xlabel('Period')
ylabel('Power')
end


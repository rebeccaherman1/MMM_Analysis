function [ft, freqs, positive_periodogram, positive_frequencies ] = powerspectrum( hi, ti, add_to_fig, use_period, spec)
%POWERSPECTRUM hi = sequence, a COLUMN vector, or the COLUMNS of a matrix. ti = time
%returns periodogram and frequencies. also plots periodogram vs.
%frequencies. Use add_to_fig=false to add to previously existing figure.
if(nargin < 3)
    add_to_fig = false;
end
if(nargin < 4)
    use_period = false;
end
if(nargin < 5)
    special_specs = false;
else
    special_specs = true;
end
N = length(ti);
ft = fft(hi); %uses columns!!!!
df = 1/(ti(end)-ti(1));
%fmax = df*N/2;
num_freqs = floor(N/2); % = fmax/df;
freqs = [0:num_freqs, -num_freqs+1:-1]'*df;
positive_frequencies = freqs(2:num_freqs)/sqrt(N);%necessary to satisfy Parseval's theorem
periodogram = ft.*conj(ft);
positive_periodogram = periodogram(2:num_freqs,:);
if(sum(positive_periodogram==0)>0)
    positive_periodogram(positive_periodogram==0) = min(positive_periodogram(positive_periodogram~=0))/10^35;
end
positive_ft = ft(2:num_freqs,:);
if(~add_to_fig)
    clf()
else
    hold on;
end
if(use_period)
    loglog(1./positive_frequencies, positive_periodogram)
    hold on
    xlabel('Period')
else
    if(special_specs)
        loglog(positive_frequencies, positive_periodogram, spec, 'Linewidth', 2);
    else
        semilogy(positive_frequencies, positive_periodogram)
    end
    xlabel('Frequency')
end
ylabel('Power')
end


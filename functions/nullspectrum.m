function [null, f] = nullspectrum(hi, periodogram, use_period, smooth_n, around_null)
%POWERSPECTRUM hi = sequence, a COLUMN vector, or the COLUMNS of a matrix. ti = time
%returns periodogram and frequencies. also plots periodogram vs.
%frequencies. Use add_to_fig=false to add to previously existing figure.
if(nargin < 3)
    use_period = false;
end
if(nargin < 4)
    smooth_n = 1;    
end
if(nargin < 5)
    around_null = true;
end
    
[null, f] = pcov(hi-mean(hi), 1, length(hi));
%The mean was removed here because pcov behaves strangely when the mean
%isn't first removed. I investigated by constructing the null spectrum by
%hand, and the true null is not affected by the presence of a mean, and
%looks like the output of pcov when the mean is first removed.
p = sum(periodogram(2:end));
null = null/sum(null)*p/2;
f = f/2/pi;
dof = 2*smooth_n;
CL_spec = dof/chi2inv(.95, dof);
CU_null = chi2inv(.95, dof)/dof;
hold on
if(use_period)
    loglog(1./f, null, 'DisplayName', 'Null Spectrum')
    xlabel('Period')
    if(around_null)
        loglog(1./f, null*CU_null, 'DisplayName', '95% confidence level')
    else
        loglog(1./f, periodogram*CL_spec, 'DisplayName', '95% confidence level')
    end
else
    semilogy(f, null, 'DisplayName', 'Null Spectrum')
    xlabel('Frequency')
    if(around_null)
        semilogy(f, null*CU_null, 'DisplayName', '95% confidence level')
    else
        semilogy(f, null*CL_spec, 'DisplayName', '95% confidence level')
    end
end
end

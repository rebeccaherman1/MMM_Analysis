NARI = load('data/ts/observations.mat');
NARI = NARI.var(:,:,3);
amip_piF = load('Analysis/pr/amip-piF_1901-2014_N500.mat');
obs = load('data/pr/observations.mat');
obs = obs.var(1901<=obs.T & obs.T <= 2014);
amip_hist = load('Analysis/pr/amip-hist_1901-2014_N500.mat');

NARI_lo = smooth(NARI, 19);

corr(NARI', obs')
corr(NARI', amip_hist.MMM.MMM')

r = corr(NARI', amip_piF.MMM.MMM')
rs = corr(NARI', amip_piF.historical_bootstrapped.b_means');
[r_lo, r_hi] = confidence_interval(rs,2);
r_hat = mean(rs)
r_mid = mean([r_lo, r_hi])
r_2sig_hat = 2*std(rs)
r_2sig_mid = mean([r_mid-r_lo, r_hi-r_mid])

C = cov(NARI', amip_piF.MMM.MMM');
t = C(1,2)/C(1,1)
ts = nan(500,1);
for i = 1:500
    bm = amip_piF.historical_bootstrapped.b_means(i,:);
    C = cov(NARI', bm');
    ts(i) = C(1,2)/C(1,1);
end
[t_lo, t_hi] = confidence_interval(ts);
t_hat = mean(ts)
t_mid = mean([t_lo, t_hi])
t_2sig_hat = 2*std(ts)
t_2sig_med = mean([t_mid - t_lo, t_hi - t_mid])

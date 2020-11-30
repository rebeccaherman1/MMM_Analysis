SPAN = 19;

AA5 = load('Analysis/pr/a_1901-2003_N500.mat');
AA6 = load('Analysis/pr/cmip6_a_1901-2014_N500.mat');
AA5_lo = smooth(AA5.MMM.MMM, SPAN)';
AA5_hi = AA5.MMM.MMM - AA5_lo;
AA6_lo = smooth(AA6.MMM.MMM, SPAN)';
AA6_hi = AA6.MMM.MMM - AA6_lo;
obs = load('data/pr/observations.mat');
obs5 = obs.var(1901<=obs.T & obs.T<=2003);
obs6 = obs.var(1901<=obs.T & obs.T<=2014);

R = corr([AA5_lo', AA5_hi', AA6_lo(1:103)', AA6_hi(1:103)', obs5']);
C = cov([AA5_lo', AA5_hi', AA6_lo(1:103)', AA6_hi(1:103)', obs5']);

coefs5 = ([C(1,1), C(2,2)]/sum([C(1,1), C(2,2), 2*C(1,2)])).^.5
coefs6 = ([C(3,3), C(4,4)]/sum([C(3,3), C(4,4), 2*C(3,4)])).^.5

R_full = corr([AA5.MMM.MMM', AA6.MMM.MMM(1:103)', obs5']);

change_lo = coefs5(1)*R(end, 1) - coefs6(1)*R(end, 3)
change_hi = coefs5(2)*R(end, 2) - coefs6(2)*R(end, 4)
change_tot = R_full(end,1)-R_full(end,2)

change_lo/change_tot
change_hi/change_tot
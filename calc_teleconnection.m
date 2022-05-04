end_yr = 2014;

obs = load('data/pr/7-9/observations.mat');
obs = obs.var(1901<=obs.T & obs.T <= end_yr);
NARI = load('data/ts/7-9/observations.mat');
T_use = NARI.T >= 1901 & NARI.T <= end_yr;
NARI = NARI.var(:,T_use,strcmp(NARI.indices, 'NARI'));

C = cov(NARI', obs');
t = C(1,2)/C(1,1);
fprintf('For obs: %f\n\n', t)

sims = {'amip-piF', 'amip-hist', 'cmip6_h', 'cmip6_a', 'cmip6_n', 'cmip6_g'};
%sims = {'h', 'a', 'n', 'g'};
t_values = nan(1, length(sims));
for s = 1:length(sims)
    sm = sims{s};
    S = load(['Analysis/pr/',sm,'_1901-',num2str(end_yr),'_N500.mat']);
    if(~contains(sm, 'amip'))
        NARI = load(['Analysis/ts/',sm,'_1901-',num2str(end_yr),'_N500.mat']);
        NARI = NARI.MMM.MMM(:,:,3);
    end
    
    r = corr(NARI', S.MMM.MMM');
    rs = corr(NARI', S.historical_bootstrapped.b_means');
    [r_lo, r_hi] = confidence_interval(rs,2);
    r_hat = mean(rs);
    r_mid = mean([r_lo, r_hi]);
    r_2sig_hat = 2*std(rs);
    r_2sig_mid = mean([r_mid-r_lo, r_hi-r_mid]);
    
    C = cov(NARI', S.MMM.MMM');
    t = C(1,2)/C(1,1);
    t_values(s) = t;
    ts = nan(500,1);
    for i = 1:500
        bm = S.historical_bootstrapped.b_means(i,:);
        C = cov(NARI', bm');
        ts(i) = C(1,2)/C(1,1);
    end
    [t_lo, t_hi] = confidence_interval(ts);
    t_hat = mean(ts);
    t_mid = mean([t_lo, t_hi]);
    t_2sig_hat = 2*std(ts);
    t_2sig_med = mean([t_mid - t_lo, t_hi - t_mid]);
    fprintf('For %s:\n %f, [%f, %f]\n %f +/- %f\n %f +/- %f\n\n', sm, t, t_lo, t_hi, t_hat, t_2sig_hat, t_mid, t_2sig_med)
end

H = load('data/pr/7-9/cmip6_a_GM.mat');
indiv_pr = load('data/pr/7-9/cmip6_piC_all.mat');
indiv_ts = load('data/ts/7-9/cmip6_piC_all.mat');
%prm = arrayfun(@(rw) strjoin(indiv_pr.model(rw,:), ' '), (1:size(indiv_pr.model,1))', 'UniformOutput', false);
%tsm = arrayfun(@(rw) strjoin(indiv_ts.model(rw,:), ' '), (1:size(indiv_ts.model,1))', 'UniformOutput', false);
PT = table(indiv_pr.model(:,1), indiv_pr.model(:,2), indiv_pr.model(:,3), indiv_pr.runs, 'VariableNames', {'institution', 'model', 'run', 'pr'}); 
TT = table(indiv_ts.model(:,1), indiv_ts.model(:,2), indiv_ts.model(:,3), indiv_ts.runs(:,:,3), 'VariableNames', {'institution', 'model', 'run', 'ts'});
TT = TT(ismember(TT.institution, H.models(:,1)),:);
T_both = innerjoin(PT, TT, 'Keys', {'institution', 'model', 'run'});
T_both.pr(T_both.pr==0)=nan;
T_both.ts(T_both.ts==0)=nan;
%what models are missing?
H.models(~ismember(H.models, T_both.model),:)
%{
pr_inc = contains(prm, tsm); [prM, prI] = sort(prm(pr_inc));
ts_inc = contains(tsm, prm); [tsM, tsI] = sort(tsm(ts_inc));
pr_t = indiv_pr.time(pr_inc,:); ts_t = indiv_ts.time(ts_inc,:); 
dimlim = min(size(pr_t,2), size(ts_t,2));
pr_t = pr_t(prI,1:dimlim); ts_t = ts_t(tsI,1:dimlim);
pr = indiv_pr.runs(pr_inc,:); pr = pr(prI,1:dimlim); 
pr(pr==0)=nan; 
ts = indiv_ts.runs(ts_inc,:,3); ts = ts(tsI,1:dimlim); 
ts(ts==0)=nan; 
%}
t_indiv = nan(height(T_both),1);
for i = 1:height(T_both)
    C = nancov(T_both.ts(i,:)',T_both.pr(i,:)');
    t_indiv(i) = C(1,2)/C(1,1);
end
T_both.tc = t_indiv;

[model_names, ~, model_groupings] = unique(T_both.model); 
MM_tc = splitapply(@mean, T_both.tc, model_groupings);
[institution_names, ~, institution_groupings] = unique(model_names);
IM_tc = splitapply(@mean, MM_tc, institution_groupings);
TC = table(institution_names, IM_tc, 'VariableNames', {'institution', 'Teleconnection'});
[~,I] = sort(abs(TC.Teleconnection-t_values(1)));
TC(I,:)

[lo, hi] = confidence_interval(t_indiv);
fprintf('For piC simulations:\n %f [%f, %f]\n %f +/- %f\n %f +/- %f\n\n', mean(t_indiv), lo, hi, mean(t_indiv), 2*std(t_indiv), mean([lo, hi]), mean([mean([lo, hi])-lo, hi-mean([lo, hi])]))
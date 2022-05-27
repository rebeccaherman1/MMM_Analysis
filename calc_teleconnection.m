%sims = {'amip-piF', 'amip-hist', 'cmip6_h', 'cmip6_a', 'cmip6_n', 'cmip6_g'}; end_yr = 2014;
sims = {'h', 'a', 'n', 'g'}; end_yr = 2003;


obs = load('data/pr/7-9/observations.mat');
obs = obs.var(1901<=obs.T & obs.T <= end_yr);
NARI = load('data/ts/7-9/observations.mat');
T_use = NARI.T >= 1901 & NARI.T <= end_yr;
NARI = NARI.var(:,T_use,strcmp(NARI.indices, 'NARI'));

C = cov(NARI', obs');
t = C(1,2)/C(1,1);
fprintf('For obs: %f\n\n', t)


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

A = load('data/pr/7-9/cmip6_a_GM.mat');
H = load('data/pr/7-9/cmip6_h_MM.mat');
indiv_pr = load('data/pr/7-9/cmip6_piC_all.mat');
indiv_ts = load('data/ts/7-9/cmip6_piC_all.mat');
indivs = {indiv_pr, indiv_ts};
for i = 1:length(indivs)
    indiv = indivs{i};
    if(any(strcmp(indiv.model(:,1), 'CESM')))
        indiv.model(strcmp(indiv.model(:,1), 'CESM'),1)={'NCAR'};
    end
    if(any(strcmp(indiv.model(:,1), 'CanESM')))
        indiv.model(strcmp(indiv.model(:,1), 'CanESM'),1)={'CCCma'};
    end
    if(any(strcmp(indiv.model(:,1), 'GISS')))
        indiv.model(strcmp(indiv.model(:,1), 'GISS'),1)={'NASA'};
    end
    if(any(strcmp(indiv.model(:,1), 'NorESM')))
        indiv.model(strcmp(indiv.model(:,1), 'NorESM'),1)={'Nor'};
    end
end


%prm = arrayfun(@(rw) strjoin(indiv_pr.model(rw,:), ' '), (1:size(indiv_pr.model,1))', 'UniformOutput', false);
%tsm = arrayfun(@(rw) strjoin(indiv_ts.model(rw,:), ' '), (1:size(indiv_ts.model,1))', 'UniformOutput', false);
PT = table(indiv_pr.model(:,1), indiv_pr.model(:,2), indiv_pr.model(:,3), indiv_pr.runs, 'VariableNames', {'institution', 'model', 'run', 'pr'}); 
TT = table(indiv_ts.model(:,1), indiv_ts.model(:,2), indiv_ts.model(:,3), indiv_ts.runs(:,:,3), 'VariableNames', {'institution', 'model', 'run', 'ts'});
TT = TT(ismember(TT.institution, A.models(:,1)),:);
TT = TT(ismember(TT.model, H.models(:,2)),:);
T_both = innerjoin(PT, TT, 'Keys', {'institution', 'model', 'run'});
T_both.pr(T_both.pr==0)=nan;
T_both.ts(T_both.ts==0)=nan;
%what models are missing?
A.models(~ismember(A.models, T_both.model),:)
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

[model_names, I, model_groupings] = unique(T_both.model); 
MM_tc = splitapply(@mean, T_both.tc, model_groupings);
[institution_names, ~, institution_groupings] = unique(T_both.institution(I));
IM_tc = splitapply(@mean, MM_tc, institution_groupings);
TC = table(institution_names, IM_tc, 'VariableNames', {'institution', 'Teleconnection'});
[~,I] = sort(abs(TC.Teleconnection-t_values(1)));
TC(I,:)

HA = load(['Analysis/pr/',sims{end-3},'_1901-',num2str(end_yr),'_N500.mat']);
stts = table(HA.indiv_s.models, HA.indiv_s.r, HA.indiv_s.e,...
    'VariableNames', {'institution', 'r', 'e'})
H_ts = load('data/ts/7-9/cmip6_h_GM.mat');
H_pr = load('data/pr/7-9/cmip6_h_GM.mat');
H_tc = nan(size(H_ts.GMs,1),1);
for i = 1:size(H_ts.GMs,1)
    ts_m = H_ts.GMs(i,:,3);
    p_m = H_pr.GMs(i,:);
    C = cov(ts_m', p_m');
    H_tc(i) = C(1,2)/C(1,1);
end
H_tc = table(H_ts.models, H_tc, 'VariableNames', {'institution', 'hist-tel'});

J = innerjoin(TC, stts);
J = join(J, H_tc)
figure()
subplot(1,2,1)
scatter(J.Teleconnection, J.r, 'filled')
xlabel('NARI Teleconnection Strength')
ylabel('r_L_F')
hold on; yl = ylim; plot(t_values(1)*[1,1], yl, 'k-');
subplot(1,2,2)
scatter(J.Teleconnection, J.e, 'filled')
xlabel('NARI Teleconnection Strength')
ylabel('sRMSE_L_F')
hold on; yl = ylim; plot(t_values(1)*[1,1], yl, 'k-');

[lo, hi] = confidence_interval(MM_tc);
fprintf('For piC simulations:\n %f [%f, %f]\n %f +/- %f\n %f +/- %f\n\n', mean(t_indiv), lo, hi, mean(t_indiv), 2*std(t_indiv), mean([lo, hi]), mean([mean([lo, hi])-lo, hi-mean([lo, hi])]))
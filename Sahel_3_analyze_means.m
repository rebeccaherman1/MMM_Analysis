clear
N = 500;
short=false;

dt = "";%, "detrended"];
%fl = "last";%, "first"];
scenarios = {'v'};%'r'};%'a6'};%'e'};%'h'};%,'a','n','g'};%'amip',; 

global start_year; 
if(short)
    start_year = 1950;% 
else
    start_year = 1901;%
end
global end_year; end_year = 2003;
global ref_T_years ref_T_amip

mkdir('analysis')

obs = load('model_output/historical_precipitation.mat');
obs_anomaly = obs.prcp-mean(obs.prcp);
ref_T_years = obs.T; ref_T_amip = start_year:end_year;

for j = 1:length(scenarios)
    scenario = scenarios{j};
    fprintf("Accessing scenario %s\n", scenario);

    h = load(['model_output/', scenario,'_GM.mat']);
    num_models = length(h.models); %num_pC_models = length(h.piC_models); 
    fname = ['analysis/', scenario, '_N', num2str(N)];
    if(short)
        Analysis = matfile([fname, '_', num2str(start_year)],'Writable', true);%'_short'], 
    else
        Analysis = matfile(fname, 'Writable', true);
    end

    T_obs = ref_T_years>=start_year & ref_T_years<=end_year;
    if(strcmp(scenario, 'amip')==2 || short)
        T_m = ref_T_amip>=start_year & ref_T_amip<=end_year;
    else
        T_m = T_obs;
    end
    hm = h.GMs(:,T_m); %TODO should actually add a T record to the data I save :/
    trust = h.trust; o = obs_anomaly(T_obs); 
   % pC = h.piC_GM; pC_trust = h.piC_trust; %(:,~isnan(sum(h.piC_GM, 1))) removes columns with NaN

    [r, e, mmm] = calc_stats(hm, trust, o);
    MMM.r = r; MMM.e = e; MMM.MMM = mmm; Analysis.MMM = mmm;
    
    hall = load(['model_output/', scenario, '_all.mat']);
    runs = hall.runs(:,T_m);
    runs_r = m_corrcoef(runs, o); runs_e = rmse(runs, o);
    indiv_runs.r = runs_r; indiv_runs.e = runs_e; indiv.models = hall.model;
    
    ir = m_corrcoef(hm, o); ie = rmse(hm, o);
    indiv.r = ir; indiv.e = ie; indiv.models = h.models;
    [~, r_order] = sort(ir, 'descend'); best_r = h.models(r_order);
    [~, e_order] = sort(ie, 'ascend'); best_e = h.models(e_order);
    indiv.best_models_r = best_r; indiv.best_models_e = best_e;
    
    [r_sanity, e_sanity, ~] = calc_stats(zeros(size(o)), 1, o);
    sanity.r = r_sanity; sanity.e = e_sanity;
    
    Analysis.MMM = MMM; Analysis.indiv = indiv; Analysis.sanity = sanity; Analysis.N = N; Analysis.indiv_runs = indiv_runs;
      
    Analysis.historical_bootstrapped = bootstrap_model(N, obs_anomaly, T_obs, hm, trust);
    %Analysis.piC_resampled_bootstrapped = sample_model(N, obs_anomaly, T_obs, pC, pC_trust, dt);
    
    %Analysis.piC_last = sample_model(N, obs_anomaly, T_obs, pC, pC_trust, dt, 'last');
end

%T can be a logical array or a range.
%TODO put detrending in here
function [r, e, mmm] = calc_stats(means, trust, obs)
    weights = trust / mean(trust);
    %mmm_std = mean(std(means,0,2));
    mmm = mean(weights.*means,1);
    e = rmse(mmm, obs);
    R = corrcoef(mmm, obs);
    r = R(1,2);
    if(isnan(r))
        r = 0;
    end
end

function [all] = bootstrap_model(N, obs, T, GMs, trust, fl, dt)
    global start_year end_year
    num_models=length(trust);
    if(nargin<6)
        means = GMs;%(:,T); TODO
    elseif(strcmp(fl, "first"))
        means = GMs(:,1:(end_year-start_year+1));
    elseif(strcmp(fl, "last"))
        length_of_GMs = sum(~isnan(a),2);
        idx_raw = cumsum(ones(size(GMs)), 2);
        idx = (idx_raw > length_of_GMs - sum(T)) & (idx_raw <= length_of_GMs);
        m_t = mod.means'; idx_t = idx';
        means = vec2mat(m_t(idx_t), sum(T));
    else
        fprintf("unrecognized option for 'fl': %s", fl);
        rs = NA; es = NA; 
        return
    end
    means = means - mean(means, 2);
    if(nargin>6 && strcmp(dt, "detrended"))
        means = detrend(means')';
    end

    rs = zeros(N,1); es = zeros(N,1);
    mmms = zeros(N, end_year-start_year+1);
    for bs = 1:N
        index=randi(num_models,num_models,1);
        reordered_means=means(index,:);
        reordered_trust=trust(index);
        [r, e, mmm] = calc_stats(reordered_means, reordered_trust, obs(T));
        mmms(bs,:) = mmm; 
        rs(bs) = r; es(bs) = e;
    end
    [CI_low, CI_high] = confidence_interval(mmms);
    b_means = mmms; mids = mean(mmms, 1);
    
    all.b_means = b_means; all.rs = rs; all.es = es; all.high = CI_high; all.low = CI_low;
    all.N = N; all.num_models = num_models; all.mids = mids;
end

function [all] = sample_model(N, obs, T, GMs, trust, dt, fl)
    num_models=length(trust);
    idx_raw = cumsum(ones(size(GMs)), 2);
    length_run = sum(T);
    length_of_GMs = sum(~isnan(GMs),2);

    rs = zeros(N,1); es = zeros(N,1); 
    mmms = zeros(N, length_run);
    for bs = 1:N
        if (nargin>6 && strcmp(fl, "last"))
            rand_end_offset = length_of_GMs;
        else
            rand_end_offset = zeros(num_models, 1);
            for mn = 1:num_models
                rand_end_offset(mn) = randi([length_run, length_of_GMs(mn)]);
            end    
        end
        idx = (idx_raw > rand_end_offset - length_run) & (idx_raw <= rand_end_offset);
        m_t = GMs'; idx_t = idx';
        means = vec2mat(m_t(idx_t), length_run);
        means = means - mean(means, 2);
        
        if(nargin>5 && strcmp(dt, "detrended"))
            means = detrend(means')';
        end
        %still bootstrap after
        index=randi(num_models,num_models,1);
        means=means(index,:);

        [r, e, mmm] = calc_stats(means, trust(index), obs(T));
        mmms(bs,:) = mmm;
        rs(bs) = r; es(bs) = e; 
    end
    [CI_low, CI_high] = confidence_interval(mmms); mids = mean(mmms, 1);
    r_means = mmms;
    %r_means, rs, es, CI_low, CI_high
    all.r_means = r_means; all.rs = rs; all.es = es; all.low = CI_low; all.high = CI_high;
    all.num_models = num_models; all.mids = mids;
end

%takes the correlation coefficient of horizontal time series in a 2D matrix 
%DATA_M with a single time series in a horizontal vector data_v. Returns a
%vertical vector of Correlation Coefficients.
function [cc] = m_corrcoef(data_m, data_v)
    n = length(data_v);
    mu_m = mean(data_m, 2); mu_v = mean(data_v);
    %"1" signifies that we will use N instead of N-1.
    sigma2_m = var(data_m,1,2); sigma2_v = var(data_v,1);
    cc = 1/n*(data_m-mu_m)*(data_v'-mu_v)./(sigma2_m.*sigma2_v).^(.5);
end

function [error] = rmse(data_m, data_v)
    data_m = data_m-mean(data_m, 2); data_v = data_v - mean(data_v);
    error = mean((data_m - data_v).^2,2).^(.5)/std(data_v,1);
end    

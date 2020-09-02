clear
N = 500;

dt = "";%, "detrended"];
%fl = "last";%, "first"];
scenarios = {'v'};%'cmip6_h', 'cmip6_a', 'cmip6_n', 'cmip6_g'};%'r'};%'a6'};%'e'};%'h'};%,'a','n','g'};%'amip',; 
variables = {'pr'};%'pr', 
realm = 'amip';

global start_year end_year ref_T_years
start_year = 1901;
end_year = 2014; 

mkdir('analysis')

for v = 1:length(variables)
    variable = variables{v};
    mkdir(['analysis/', variable])
    obs = load(['data/', variable, '/observations.mat']);
    obs_anomaly = obs.var-mean(obs.var);
    timeframe_obs = (obs.T >= start_year & obs.T <= end_year);
    ref_T_years = obs.T(timeframe_obs);
    o = obs_anomaly(:,timeframe_obs,:);
    for j = 1:length(scenarios)
        scenario = scenarios{j};
        fprintf("Accessing scenario %s\n", scenario);

        h = load(['data/', variable, '/', scenario,'_GM.mat']);
        fname = ['analysis/', variable, '/', scenario, '_', num2str(ref_T_years(1)), '-', num2str(ref_T_years(end)), '_N', num2str(N)];
        hall = load(['data/', variable, '/', scenario, '_all.mat']);
        num_models = length(h.models); %num_pC_models = length(h.piC_models); 
        Analysis = matfile(fname, 'Writable', true);

        T_m = ismember(single(h.time(1,:)),ref_T_years);
        hm = h.GMs(:,T_m, :); %TODO should actually add a T record to the data I save :/
        trust = h.trust;

        [r, e, mmm] = calc_stats(hm, trust, o);
        MMM.r = r; MMM.e = e; MMM.MMM = mmm; Analysis.MMM = mmm;

        %hall
        runs = hall.runs(:,T_m,:);
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

        Analysis.historical_bootstrapped = bootstrap_model(N, o, timeframe_obs, hm, trust);

        %TODO could alternately use ISFIELD and then I wouldn't have to
        %define the realm at the top...
        if(~strcmp(realm, 'amip'))
            [Analysis.piC_resampled_bootstrapped, skip_models] = sample_model(N, o, timeframe_obs, h.piC_GMs, h.piC_trust, dt);
            fprintf('skipping piC simulations which are too short:')
            h.piC_models(skip_models,:)
        end
        
        if(isfield(h, 'indices'))
            Analysis.indices = h.indices;
        end
    end
end

%% functions

%T can be a logical array or a range.
%TODO put detrending in here
function [r, e, mmm] = calc_stats(means, trust, obs)
    weights = trust / sum(trust);
    %mmm_std = mean(std(means,0,2));
    mmm = sum(weights.*means,1);
    e = rmse(mmm, obs);
    [~,~,s3]=size(means);
    r=nan([1,1,s3]);
    for i=1:s3
        r(:,:,i) = corr(mmm(:,:,i)', obs(:,:,i)');
    end
end

function [all] = bootstrap_model(N, obs, T, GMs, trust, fl, dt)
    run_length = sum(T);
    num_models=length(trust);
    [s1,s2,s3]=size(GMs);
    if(nargin<6)
        means = GMs;%(:,T); TODO
    elseif(strcmp(fl, "first"))
        means = GMs(:,1:(end_year-start_year+1),:);
    elseif(strcmp(fl, "last"))
        length_of_GMs = sum(~isnan(GMs(:,:,1)),2);
        idx_raw = cumsum(ones([s1,s2]), 2);
        idx = (idx_raw > length_of_GMs - run_length) & (idx_raw <= length_of_GMs);
        m_t = mod.means'; idx_t = idx';
        %technically I should put this in a for loop, but I don't really
        %care!!!!! because I don't use first and last anymore.
        means = vec2mat(m_t(idx_t), run_length);
    else
        fprintf("unrecognized option for 'fl': %s", fl);
        rs = NA; es = NA; 
        return
    end
    means = means - mean(means, 2);
    if(nargin>6 && strcmp(dt, "detrended"))
        means = permute(detrend(permute(means,[2,1,3])), [2,1,3]);
    end

    rs = zeros(N,1,s3); es = zeros(N,1,s3);
    mmms = zeros(N, run_length,s3);
    for bs = 1:N
        index=randi(num_models,num_models,1);
        reordered_means=means(index,:,:);
        reordered_trust=trust(index);
        [r, e, mmm] = calc_stats(reordered_means, reordered_trust, obs(:,T,:));
        mmms(bs,:,:) = mmm; 
        rs(bs,:,:) = r; es(bs,:,:) = e;
    end
    [CI_low, CI_high] = confidence_interval(mmms);
    b_means = mmms; mids = mean(mmms, 1);
    
    all.b_means = b_means; all.rs = rs; all.es = es; all.high = CI_high; all.low = CI_low;
    all.N = N; all.num_models = num_models; all.mids = mids;
end

function [all, too_short] = sample_model(N, obs, T, GMs, trust, dt, fl)
    length_run = sum(T);
    length_of_GMs = sum(~isnan(GMs(:,:,1)),2);
    too_short = length_of_GMs<length_run;
    
    length_of_GMs = length_of_GMs(~too_short);
    GMs = GMs(~too_short,:,:); trust = trust(~too_short);
    [s1,s2,s3] = size(GMs);
    idx_raw = cumsum(ones([s1,s2]), 2);
    
    num_models=length(trust);

    rs = zeros(N,1,s3); es = zeros(N,1,s3); 
    mmms = zeros(N, length_run,s3);
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
        m_t = permute(GMs, [2,1,3]); idx_t = repmat(idx',1,1,s3);
        means = permute(reshape(m_t(idx_t), length_run,s1,s3), [2,1,3]);
        means = means - mean(means, 2);
        
        if(nargin>5 && strcmp(dt, "detrended"))
            %this won't work with SSTs!
            means = detrend(means')';
        end
        %still bootstrap after
        index=randi(num_models,num_models,1);
        means=means(index,:,:);

        [r, e, mmm] = calc_stats(means, trust(index), obs(:,T,:));
        mmms(bs,:,:) = mmm;
        rs(bs,:,:) = r; es(bs,:,:) = e; 
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
    [~,n,~] = size(data_v);
    mu_m = mean(data_m, 2); mu_v = mean(data_v,2);
    %"1" signifies that we will use N instead of N-1.
    sigma2_m = var(data_m,1,2); sigma2_v = var(data_v,1,2);
    cc = 1/n*sum((data_m-mu_m).*(data_v-mu_v),2)./(sigma2_m.*sigma2_v).^(.5);
end

function [error] = rmse(data_m, data_v)
    data_m = data_m-mean(data_m, 2); data_v = data_v - mean(data_v,2);
    error = mean((data_m - data_v).^2,2).^(.5)./std(data_v,1,2);
end    

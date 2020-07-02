clear; 
scenarios = {'v'};%'r'};%'a6'};%'e';%'hgna'; 
single_scenario = scenarios{1};
colors = 'bgrm'; 
names = [{'v'}];%'r'}];%'a6'}];%[{'e'}];%[{'ALL'}, {'GHG'}, {'NAT'}, {'AA'}];
global end_year start_year PAD; 
start_year = 1901;%50; 
end_year = 2003; 
PAD = 1101;

wetdry = [   84,48,5;  140,81,10;191,129,45;223,194,125;218,204,167;189,189,189
         171,206,201;128,205,193;53,151,143;   1,102,94;    0,60,48]/255;
x = round(linspace(0, 255, 11))';
map = interp1(x/255,wetdry,linspace(0,1,255));
clear wetdry x;

obs = load('data/historical_precipitation.mat'); 
obsT = obs.T(obs.T>=start_year & obs.T<=end_year)'; obsP_o = obs.prcp(obs.T>=start_year & obs.T<=end_year)'; clear obs;
obsP = [obsP_o;mean(obsP_o)*ones(PAD-length(obsP_o),1)];


%MAKE FIGURE 3
%TODO can make similar splitting change here
%also perhaps can make match number of runs in included models exactly (For
%an under- instead of over-estimate).
figure(3); clf; set(gca, 'fontsize', 15)

for i=1:length(scenarios)
    s = scenarios{i};
    c = colors(i);%scenarios==s);
    aname = ['analysis/', s, '_N500'];
    if(start_year>1901)
	aname=[aname, '_', num2str(start_year)];
    end
    A = load([aname, '.mat']);
    forced_mmms = A.historical_bootstrapped.b_means; sz = size(forced_mmms);
    %the analysis file actually already has the mean subtracted. Here was the bug.
    forced_mmm_anomalies = (forced_mmms - mean(forced_mmms,2));
    [periods, mp, low, high, ~, ~] = ps_multi(forced_mmm_anomalies', obsT, .17, PAD);
    fill([periods;flipud(periods)],[low';flipud(high')],c,'FaceAlpha', .3 ,'linestyle','none', 'DisplayName', ['95% CI ', names{i}, ' PS']);%scenarios==s
    hold on;
    %{
    if(s=='a')
        unforced_mmms = A.piC_resampled_bootstrapped.r_means;
        unforced_mmm_anomalies = unforced_mmms - mean(unforced_mmms,2);
        [periods, ~, low, high, ~, ~] = ps_multi(unforced_mmm_anomalies, obsT, .17,PAD);
        fill([periods;flipud(periods)],[low';flipud(high')],'y','FaceAlpha', .3 ,'linestyle','none', 'DisplayName', '95% CI piC PS');
    end
    %}
    %hold on; %p.HandleVisibility = 'off';
end

for i=1:length(scenarios)
    s = scenarios{i};
    c = colors(i);%scenarios==s);
    aname = ['analysis/', s,'_N500'];
    if(start_year>1901)
        aname = [aname, '_', num2str(start_year)];
    end	
    A = load([aname, '.mat']);
    [~, ~, mean_periodogram, positive_frequencies ] = ps_mean_sahel(A.MMM.MMM', obsT, ['Mean ', names{i}, ' PS'], c, PAD);
    %{
    forced_mmms = A.historical_bootstrapped.b_means';
    forced_mmm_anomalies = (forced_mmms - mean(forced_mmms,2));
    [periods, mean_periodogram, ~, ~, ~, ~] = ps_multi(forced_mmm_anomalies, obsT, .17,PAD);
    %}
%    plot(1./positive_frequencies, mean_periodogram, [c,'-'], 'DisplayName', ['Mean ', names{i}, ' PS']);%scenarios==s
    hold on; 
    if(s=='a')
        unforced_mmms = A.piC_resampled_bootstrapped.r_means';
        unforced_mmm_anomalies = unforced_mmms - mean(unforced_mmms,1);
        [periods, mean_periodogram, ~, ~, ~, ~] = ps_multi(unforced_mmm_anomalies, obsT, .17,PAD);
        plot(periods, mean_periodogram, 'y-', 'DisplayName', 'Mean piC PS');
    end
end
legend('location', 'northwest'); xlim([2,99]);
ylabel('Power, ((mm/day)^2/frequency)');
ax = gca; ax.XScale = 'log'; %ylim([0,.3])
xlabel('Period (Years)'); title('Power Spectra of MMMs')
savefig([single_scenario, '_Fig3.fig']);

%MAKE FIGURE 5
aname = ['analysis/', single_scenario, '_N500']; %TODO! h vs e
if(start_year>1901)
    aname = [aname,'_',num2str(start_year)];end
A = load([aname,'.mat']);
mmm = A.MMM.MMM'; 
SCALE=1;%2.8;

h_all = load(['data/', single_scenario, '_all.mat']); %TODO!  h vs e
%TODO scaling
h_runs = h_all.runs'+(SCALE-1)*(mmm-mean(mmm));%(:, 1:99)'; %so individual runs will be vertical for PS

total_rainfall = mean(h_runs, 1); %why did I switch these? TODO
scaling.indices = 1:length(total_rainfall); %just a number for each element for use in indexing later
[scaling.models,idx,scaling.groupings] = unique(h_all.model(:,2), 'stable'); %list of unique models; groupings for each model
[~,~,scaling.groupings2] = unique(h_all.model(idx,1)); 
scaling.models = scaling.models'; scaling.groupings = scaling.groupings'; %so that models are listed horizontally
scaling.groupings2 = scaling.groupings2';
scaling.mMTR = splitapply(@mean, total_rainfall, scaling.groupings);
half = max(abs(scaling.mMTR-mean(obsP)));
scaling.color_index_mean = floor(1+254/2*(1+(scaling.mMTR-mean(obsP))/half));
scaling.basic_name = 'a. Observations and All simulations';
normalized_runs = (h_runs)./scaling.mMTR(scaling.groupings)*mean(obsP);
normalized_runs_anomalies = normalized_runs - mean(normalized_runs,1);
clear half total_rainfall;

figure(5); clf;
subplot(2,2,1)
yl = plot_all(h_runs, map, scaling, obsP, obsT); set(gca, 'FontSize', 12); 

subplot(2,2,3);
[~, ~, obs_ps, obs_f] = ps_sahel(obsP, obsT, 'Observed Spectrum', '-', 'k', 2);
%TODO ps_multi doesn't group by model! Each run counts separately. Do I want that?
[piC_periods, h_mean_periodogram, low, high, ml, mh, ~] = ps_multi(normalized_runs_anomalies, obsT, .17, PAD, scaling.groupings, scaling.groupings2); 
fill([piC_periods;flipud(piC_periods)],[low';flipud(high')],[.8, .8, 1],'FaceAlpha', .3 ,'linestyle','none', 'DisplayName', '95% range ALL spectra by model');
hold on; set(gca, 'FontSize', 12);
fill([piC_periods;flipud(piC_periods)],[ml';flipud(mh')],[.7, .7, .9],'FaceAlpha', .3 ,'linestyle','none', 'DisplayName', '66% range ALL spectra by model');
hold on; 
for i=1:length(scenarios)
    s = scenarios{i};
    c = colors(i);%scenarios==s);
    A = load(['data/', s, '_all.mat']);
    %TODO scaling
    [~,l_idx,local_groupings] = unique(A.model(:,2), 'stable');
    [~,~,local_groupings2] = unique(A.model(l_idx,1), 'stable');
    forced_runs = A.runs'+(SCALE-1)*(mmm-mean(mmm));%(:,1:99)';
    normalized_forced_runs = forced_runs./mean(forced_runs,1)*mean(obsP); %I'm not grouping by model for this mean.
    normalized_forced_runs_anomalies = normalized_forced_runs - mean(normalized_forced_runs,1);
    [piC_periods, mean_periodogram, ~, ~, ~, ~] = ps_multi(normalized_forced_runs_anomalies, obsT, .17, PAD, local_groupings', local_groupings2');
    %f = fill([periods;flipud(periods)],[low';flipud(high')],c,'FaceAlpha', .3 ,'linestyle','none', 'DisplayName', '95% confidence Interval');
    hold on;  
    plot(piC_periods, mean_periodogram, [c,'--'], 'DisplayName', ['Tiered Mean of ', names{i}, ' spectra']);%scenarios==s
end
legend('location', 'northwest')
ylabel('Power, Scaled', 'FontWeight', 'bold'); ylim(yl); 
xlabel('Period (years)'); xlim([2,end_year-start_year+1]);
title('c.')
clear s c A forced_runs normalized_forced_runs periods mean_periodogram ans high low mh ml

clear names scenarios colors;

clear h_all h_runs

savefig([single_scenario, '_Fig5.fig']); %TODO!

%piC!
%TODO -- think about smoothing or cutting into smaller amounts and
%averaging. Does having a longer series actually not improve the error at
%all? What if I smoothed it until I had the same resolution as the
%observations?

%{
A = load('analysis/h_N500.mat');
% TODO scaling
mmm = A.MMM.MMM'; %forced = mmm - mean(mmm);
residual = obsP_o - SCALE*mmm; %forced;
residual_pad = [residual;ones(PAD-length(residual),1)*mean(residual)];
mmm_pad=[mmm;ones(PAD-length(mmm),1)*mean(mmm)];
clear A mmm

piC_all = load('data/piC_all.mat');
piC_models_old = piC_all.model; 
[to_use, ~] = ismember(piC_models_old(:,2), scaling.models); %get rid of models that don't have ALL runs
piC_runs_old = piC_all.runs(to_use,:)'; piC_models_old = piC_models_old(to_use, :);
piC_runs_old(piC_runs_old==0)=nan;
clear to_use

[piC_models, piC_runs] = split_piC_runs(piC_models_old, piC_runs_old, 103);
clear piC_runs_old normalied_piC_runs

scaling2.indices = 1:length(piC_models(:,2)); %just a number for each element for use in indexing later
[scaling2.models,idx2,scaling2.groupings] = unique(piC_models(:,2), 'stable'); %list of unique models; groupings for each model
[~,~,scaling2.groupings2] = unique(piC_models(idx2,1),'stable');
scaling2.models = scaling2.models'; scaling2.groupings = scaling2.groupings'; 
scaling2.groupings2 = scaling2.groupings2'; %so that models are listed horizontally
[~, ind] = ismember(scaling2.models, scaling.models); %match between unique piC and ALL;
scaling2.mMTR    = scaling.mMTR  (ind);
scaling2.color_index_mean = scaling.color_index_mean(ind); %floor(1 + (mMTR - min(mMTR))/rang*254);
scaling2.basic_name = 'b. Residual and piC simulations';
normalized_piC_runs = piC_runs./scaling2.mMTR(scaling2.groupings)*mean(obsP);
normalized_piC_runs_A = normalized_piC_runs - nanmean(normalized_piC_runs,1);
clear ind
subplot(2,2,4)
[piC_periods, piC_mean_periodogram, low, high, ml, mh, ~] = ps_multi(normalized_piC_runs_A, obsT, .17, PAD, scaling2.groupings, scaling2.groupings2);
[~, ~, resid_ps, ~] = ps_sahel(residual_pad, obsT, 'Residual Spectrum', '-.', 'k', 2);
fill([piC_periods;flipud(piC_periods)],[low';flipud(high')],[1, 1, .7],'FaceAlpha', .5 ,'linestyle','none', 'DisplayName', '95% range piC spectra by model');
hold on; set(gca, 'FontSize', 12);%p.HandleVisibility = 'off'; 
fill([piC_periods;flipud(piC_periods)],[ml';flipud(mh')],[.9, .9, .6],'FaceAlpha', .5 ,'linestyle','none', 'DisplayName', '66% range piC spectra by model');
%fill([periods;flipud(periods)],[(mean_periodogram-2*SD);flipud((mean_periodogram+2*SD))],[.9, .9, .7],'FaceAlpha', .3 ,'linestyle','none', 'DisplayName', '95% Confidence Interval of mean of piC spectra');%66% range ALL spectra');
hold on; %p.HandleVisibility = 'off'; 
plot(piC_periods, piC_mean_periodogram, '--', 'Color', [1, .75, 0], 'DisplayName', 'Tiered Mean of piC spectra')
legend('location', 'northwest')
xlabel('Period (years)'); xlim([2,end_year-start_year+1]); ylim(yl);
title('d.')

subplot(2,2,2)
plot_all(piC_runs, map, scaling2, residual_pad, obsT, '-.'); set(gca, 'FontSize', 12); ylim(yl);
figure(6); clf;
[~,~,mmm_periodogram,f] = ps_sahel(mmm_pad, obsT, '', '-.', 'k', 2);
piC_and_MMM = piC_mean_periodogram+mmm_periodogram;
figure(5); subplot(2,2,3);
semilogx(1./f, piC_and_MMM, 'k--', 'DisplayName', 'piC PS + ALL MMM PS');
semilogx(1./f, piC_mean_periodogram+2.8*mmm_periodogram, 'k:', 'DisplayName', 'piC PS + 2.8*ALL PS')
figure(6); clf;
plot(1./f, h_mean_periodogram, 'b--', 'DisplayName', 'Tiered Mean of ALL Spectra');
hold on;
plot(1./f, piC_and_MMM, 'k--', 'DisplayName', 'piC PS + ALL MMM PS');
legend('Location', 'northwest'); xlabel('Period (years)'); xlim([2,end_year-start_year+1]); ylabel('Power, Normalized');


%MAKE FIGURE S1
figure(6); clf; subplot(1,2,1);
plot_all(normalized_runs,     map, scaling,  obsP,     obsT, '-'); xlim([2,end_year-start_year+1]); ylim(yl);
subplot(1,2,2);
plot_all(normalized_piC_runs, map, scaling2, residual_pad, obsT, '-.'); xlim([2,end_year-start_year+1]); ylim(yl);
%clear piC_al piC_models_old piC_runs_old
clear obsP obsT residual scaling yl ans map normalized_piC_runs normalized_runs

obs_to_all = obs_ps./h_mean_periodogram; obs_all = sqrt(obs_ps) - sqrt(h_mean_periodogram);
resid_to_piC = resid_ps./piC_mean_periodogram; resid_piC = sqrt(resid_ps) - sqrt(piC_mean_periodogram);
resid_ps_interp = interp1(1./obs_f, resid_ps, piC_periods); 

piC_ps_interp = interp1(piC_periods, piC_mean_periodogram, 1./obs_f);
resid_to_piC_p = resid_ps_interp./piC_mean_periodogram; resid_piC_p = sqrt(resid_ps_interp) - sqrt(piC_mean_periodogram);
resid_to_piC_o = resid_ps./piC_ps_interp; resid_piC_o = sqrt(resid_ps) - sqrt(piC_ps_interp);
periods_op = [1./obs_f; piC_periods]; [rtp_periods, i] = sort(periods_op);
resid_to_piC_op = [resid_to_piC_o; resid_to_piC_p]; resid_piC_op = [resid_piC_o; resid_piC_p];
resid_to_piC = resid_to_piC_op(i); resid_piC = resid_piC_op(i);
clear resid_ps_interp piC_ps_interp resid_to_piC_p resid_to_piC_o periods_op resid_to_piC_op i


%MAKE FIGURE S2
figure(7); clf; 
subplot(2,1,1)
plot(1./obs_f, obs_to_all, 'b-.', 'DisplayName', 'Observed PS / mean ALL PS'); hold on; 
plot(rtp_periods, resid_to_piC, 'y-.', 'DisplayName', 'Residual PS / mean piC PS');
legend('Location', 'northwest'); xlabel('Period (years)'); xlim([2,99]); ylabel('Ratio of Normalized Power')
subplot(2,1,2)

semilogx(1./obs_f, obs_all, 'b--', 'DisplayName', 'sqrt(Observed PS) - sqrt(mean ALL PS)'); hold on; 
semilogx(1./obs_f, resid_piC, 'y--', 'DisplayName', 'sqrt(Residual PS) - sqrt(mean piC PS)');
legend('Location', 'northwest'); xlabel('Period (years)'); xlim([2,end_year-start_year+1]); ylabel('Power, Normalized');
%}
function [yl] = plot_all(m_runs, map, scaling, obsv, obsT, obs_style)
    if nargin < 6
        obs_style = '-';
    end
    global end_year start_year PAD
    %I should smooth these as well to make up for the smaller means
    for i = unique(scaling.groupings)
        ps_mean_sahel(m_runs(:, scaling.indices(scaling.groupings==i)), obsT, char(scaling.models(i)), map(scaling.color_index_mean(i),:), PAD, 1);
    end
    ps_sahel(obsv, obsT, 'Observations', obs_style, 'k', 2); %[~, ~, obs_ps, ~] = 
    %legend('show', 'Location', 'eastoutside')
    ylabel('Power, True Values', 'FontWeight', 'bold'); yl = ylim;
    xlabel(''); xlim([2,end_year-start_year+1]);
    title(scaling.basic_name, 'FontWeight', 'bold'); 
end

function [new_models, new_runs] = split_piC_runs(piC_models, piC_runs, LENGTH)
    SR = sum(floor(sum(~isnan(piC_runs),1)/LENGTH));
    new_runs = zeros(LENGTH,SR);
    new_models = cell(SR,3);
    it_max = floor(length(piC_runs(:,1))/LENGTH);
    st = 1;
    for i = 1:it_max
        columns = prod(~isnan(piC_runs((1:LENGTH)+LENGTH*(i-1),:)),1);
        new_runs(:, st:(st+sum(columns)-1)) = piC_runs((1:LENGTH)+LENGTH*(i-1),logical(columns));
        new_models(st:(st+sum(columns)-1),:) = piC_models(logical(columns'),:);
        st = st+sum(columns);
    end
end

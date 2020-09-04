clear; 

%TODO I really could use a file keeping track of the DIFFERENCE: the fast
%component.

tosave = false;
realm = 'cmip6';
variable = 'pr';
global end_year start_year PAD; 
start_year = 1901;
PAD = 1101;

switch realm
    case 'cmip6'
        scenarios = {'cmip6_h', 'cmip6_a', 'cmip6_n', 'cmip6_g'};%, 'cmip6_r', 'v'}; %add the amip ones too???
        colors = {'b', 'm', 'r', 'g', [0, 127, 0]/255, [1,.7,0]};
        names = {'CMIP6 ALL', 'CMIP6 AA', 'CMIP6 NAT', 'CMIP6 GH6', 'amip-hist', 'amip-piF'};
        end_year = 2013;
        %I COULD automate adding cmip5, but I don't feel like it...
        if strcmp(variable, 'ts')
            scenarios = scenarios(1:4); %the others will get automatically shortened. 
            end_year = 2014;
        end
    case 'cmip5'
        scenarios = {'h', 'a', 'n', 'g'};
        colors = {'b', 'm', 'r', 'g'};
        names = {'CMIP5 ALL', 'CMIP5 AA', 'CMIP5 NAT', 'CMIP5 GHG'};
        end_year = 2003;
    case 'amip'
        scenarios = {'amip-hist','amip-piF', 'cmip6_fast'};
        colors = {[0, 127, 0]/255,[1,.7,0],[126, 47, 142]/255};
        names = {'amip-hist','amip-piF', 'Implied Fast Component'};
        end_year = 2013;
    otherwise
        %a6? e? p? amip? (that last one's cmip5)
end
single_scenario = scenarios{1};

wetdry = [   84,48,5;  140,81,10;191,129,45;223,194,125;218,204,167;189,189,189
         171,206,201;128,205,193;53,151,143;   1,102,94;    0,60,48]/255;
x = round(linspace(0, 255, 11))';
map = interp1(x/255,wetdry,linspace(0,1,255));
clear wetdry x;

%horizontal vectors
obs = load(['data/', variable, '/observations.mat']); 
tm = obs.T>=start_year & obs.T<=end_year;
obsT = obs.T(tm); obsP_o = obs.var(tm); clear obs;
obsP = [obsP_o,mean(obsP_o,2)*ones(1, PAD-length(obsP_o))];
start_year = obsT(1); end_year = obsT(end);


%% MAKE FIGURE 3
%TODO can make similar splitting change here
%also perhaps can make match number of runs in included models exactly (For
%an under- instead of over-estimate).
figure(3); clf; set(gca, 'fontsize', 15)

mean_periodograms = cell(length(scenarios)+1, 1);

for i=1:length(scenarios)
    s = scenarios{i};
    c = colors{i};%scenarios==s);
    A = load(['analysis/', variable, '/', s, '_', num2str(start_year), '-', num2str(end_year), '_N500', '.mat']);
    %forced_mmms = A.historical_bootstrapped.b_means;
    
    %TO COMBINE THESE 2 SOURCES OF UNCERTAINTY:
    %go to log space: CI_log = log(CI); ps_log = log(ps); BCI = log(BCI);
    %find the std: std_log = CI_log - ps_log; B_std_log = BCI_log - ps_log;
    %convert to var: v_log = std_log.^2;
    %TODO HOW TO PRESEVE SKEWNESS?!?!
    
    %this is the version just using what pmtm gives.
    [ps, f, CI] = pmtm(A.MMM.MMM-mean(A.MMM.MMM),[],PAD, 1, 'ConfidenceLevel', .95);
    periods = 1./f(2:end); 
    
    %this is the version measuring uncertainty in the MMM via
    %bootstrapping.
    %[periods, ps, CI, skw] = pmtm_multi(A.historical_bootstrapped.b_means, obsT, [], PAD); %.17 CI?!
    %periods = periods(2:end);
    
    CI = CI(2:end,:); ps = ps(2:end); %get rid of 0/infinity
    fill([periods;flip(periods)],[CI(:,1);flip(CI(:,2))],c,'FaceAlpha', .3 ,'linestyle','none', 'HandleVisibility', 'off');%'DisplayName', ['95% CI ', names{i}, ' PS']);%scenarios==s
    mean_periodograms{i} = ps;
    hold on;
    
    if strcmp(s, 'a') || strcmp(s, 'cmip6_a')
        unforced_mmms = A.piC_resampled_bootstrapped.r_means;
        [ps, f, CI] = pmtm(unforced_mmms',[],PAD, 1, 'ConfidenceLevel', .95);
        periods = 1./f(2:end); low = mean(CI(2:end, 1:2:end), 2); high = mean(CI(2:end, 2:2:end), 2);
        %[periods, mp, low, high, ~, ~] = ps_multi(unforced_mmms, obsT, .17,PAD, 2);
        fill([periods;flip(periods)],[low;flip(high)],'y','FaceAlpha', .3 ,'linestyle','none', 'DisplayName', 'bootstrapped piC MMM PS');
        mean_periodograms{end} = mean(ps, 2);
    end
    
    %hold on; %p.HandleVisibility = 'off';
end

for i=1:length(scenarios)
    s = scenarios{i};
    c = colors{i};
    plot(periods, mean_periodograms{i}, '-', 'Color', c, 'DisplayName', [names{i}, ' MMM PS']);%scenarios==s
    hold on; 
    if strcmp(s, 'a') || strcmp(s, 'cmip6_a')
        %plot(periods, mean_periodograms{end}, 'y-', 'DisplayName', 'Mean piC PS');
    end
end
legend('location', 'northwest'); xlim([2,end_year-start_year+1]);
ylabel('Power, ((mm/day)^2/frequency)');
ax = gca; ax.XScale = 'log'; 
xlabel('Period (Years)'); title('Power Spectra of MMMs')
set(gca, 'XScale', 'linear'); xticks([2, 10:10:114])
switch realm
    case 'cmip6'
        xlim([2,50]); ylim([0,.3])
    case 'cmip5'
        xlim([2,40]); ylim([0,.2])
    case 'amip'
end
if(tosave)
    savefig(['figures/', variable, '/', realm, '_Fig3.fig']);
    saveas(3, ['figures/', variable, '/', realm, '_Fig3'], 'png');
end
%xticks([2,5,10, 20, 50, 100]); ylim([0,.3])
%(automate the tickmarks and switching to a log y-axis)
%TODO MADE it to HERE fixing my code!
%% MAKE FIGURE 5
aname = ['analysis/', variable, '/', s, '_', num2str(start_year), '-', num2str(end_year), '_N500', '.mat']; 
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

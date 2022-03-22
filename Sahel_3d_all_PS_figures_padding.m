clear; 

%TODO I really could use a file keeping track of the DIFFERENCE: the fast
%component.
%TODO only use common models...

tosave = false;
realm = 'cmip6';
variable = 'pr';
global end_year start_year PAD; 
start_year = 1901;
PAD = 1101;
start_month = 7;
end_month = 9;

isKelvin = false;
switch realm
    case 'cmip6'
        scenarios = {'cmip6_h', 'cmip6_a', 'cmip6_n', 'cmip6_g'};%, 'cmip6_r', 'v'}; %add the amip ones too???
        colors = {'b', 'm', [0.60,0.20,0.00], [0.00,0.80,0.00],[126, 47, 142]/255, [0, 127, 0]/255, [1,.7,0]};
        names = {'CMIP6 ALL', 'CMIP6 AA', 'CMIP6 NAT', 'CMIP6 GH6', 'Implied Fast Component', 'amip-hist', 'amip-piF'};
        end_year = 2013;
        %I COULD automate adding cmip5, but I don't feel like it...
        s_g = scenarios{4};
        piC_s = 'cmip6_piC';
        if strcmp(variable, 'ts')
            end_year = 2014;
            isKelvin = true;
            scenarios = scenarios(1:end-1);
        end
    case 'cmip5'
        scenarios = {'h', 'a', 'n', 'g'};
        colors = {[0.00,0.45,0.74], [0.75,0.00,0.75], 'r', [0.47,0.67,0.19]};
        names = {'CMIP5 ALL', 'CMIP5 AA', 'CMIP5 NAT', 'CMIP5 GHG'};
        end_year = 2003;
        s_g = scenarios{4};
        piC_s = 'piC';
    case 'amip'
        scenarios = {'amip-piF', 'amip-hist','cmip6_fast', 'obs'};
        colors = {[1,.7,0],[0, 127, 0]/255,[126, 47, 142]/255, 'k'};
        names = {'amip-piF', 'amip-hist','Implied Fast Component', 'Observations'};
    otherwise
        %a6? e? p? amip? (that last one's cmip5)
end
single_scenario = scenarios{1};
if(strcmp(variable, 'ts'))
    to_subtract_GHG = [false, true];
    G = load(make_analysis_filename(variable, s_g, start_year, end_year, '500'));
else
    to_subtract_GHG = [false];
end

drywet = [   84,48,5;  140,81,10;191,129,45;223,194,125;218,204,167;189,189,189
         171,206,201;128,205,193;53,151,143;   1,102,94;    0,60,48]/255;
coldhot =  [  5, 48, 97; 33,102,172; 67,147,195; 125,178,222; 183,209,240
            220,220,220; 244,193,175; 244,145,126; 214, 96, 77; 178, 24, 43
            103,  0, 31]/255;
x = round(linspace(0, 255, 11))';
if(strcmp(variable, 'pr'))
    map = interp1(x/255,drywet,linspace(0,1,255));
else
    map = interp1(x/255,coldhot,linspace(0,1,255));
end
clear wetdry hotcold x;

if(isKelvin)
    conv = permute(strcmp(variable, 'ts').*(~strcmp(G.indices, 'NARI'))*273.15,...
                   [1,3,2]);
else
    conv = 0;
end

obs = load(make_data_filename(variable, start_month, end_month, 'observations')); 
tm = obs.T>=start_year & obs.T<=end_year;
obsT = obs.T(tm); obs_var = obs.var(:,tm,:)+conv; %convert obs to Kelvin if data is in Kelvins
obs_var_a = obs_var - mean(obs_var,2);
clear obs s_g

start_year = obsT(1); end_year = obsT(end);

%% MAKE FIGURE 3
%TODO can make similar splitting change here
%also perhaps can make match number of runs in included models exactly (For
%an under- instead of over-estimate).

%for TS, I want to remove GHG from observations and other variables. 
if(strcmp(variable, 'ts'))
    I = 3;
else
    I=1;
end


figure(3); clf; set(gca, 'fontsize', 15)

for subtract_GHG = to_subtract_GHG
    mean_periodograms = cell(length(scenarios)+1, I);
    subplot(length(to_subtract_GHG),1, 1+subtract_GHG*1)

    for i=1:length(scenarios)
        s = scenarios{i};
        c = colors{i};%scenarios==s);
        if(strcmp(s, 'obs'))
            A.MMM.MMM = obs_var_a;
        elseif(subtract_GHG && contains(s, 'g'))
            continue
        else
            A = load(make_analysis_filename(variable, s, start_year, end_year, '500'));
        end
        %forced_mmms = A.historical_bootstrapped.b_means;

        I = size(A.MMM.MMM, 3); %can run through different basins if take away picking just NARI in obs above.
        for indx = I
            %subplot(1, I, indx) 

            %TO COMBINE THESE 2 SOURCES OF UNCERTAINTY:
            %go to log space: CI_log = log(CI); ps_log = log(ps); BCI = log(BCI);
            %find the std: std_log = CI_log - ps_log; B_std_log = BCI_log - ps_log;
            %convert to var: v_log = std_log.^2;
            %TODO HOW TO PRESEVE SKEWNESS?!?!

            %this is the version just using what pmtm gives.
            if(strcmp(variable, 'ts') && subtract_GHG && (contains(s, 'h') || strcmp(s, 'obs')))
                mmm = A.MMM.MMM(:,:,indx) - G.MMM.MMM(:,:,indx);
                names(i) = {[names{i}, '-GHG']};
            else
                mmm = A.MMM.MMM(:,:,indx);
            end
            [ps, f, CI] = pmtm((mmm-mean(mmm,2)),[],PAD, 1, 'ConfidenceLevel', .95);
            periods = 1./f(2:end); 

            %this is the version measuring uncertainty in the MMM via
            %bootstrapping.
            %[periods, ps, CI, skw] = pmtm_multi(A.historical_bootstrapped.b_means, obsT, [], PAD); %.17 CI?!
            %periods = periods(2:end);

            CI = CI(2:end,:); ps = ps(2:end); %get rid of 0/infinity
            if(~strcmp(realm, 'amip') && strcmp(s, 'cmip6_fast'))
                fa = 0;
                ls = '-.';
            else
                fa = .05;
                ls = ':';
            end
            fill([periods;flip(periods)],[CI(:,1);flip(CI(:,2))],c,'FaceAlpha', fa ,'linestyle',ls, 'EdgeColor', c, 'EdgeAlpha', .75,'HandleVisibility', 'off');%'DisplayName', ['95% CI ', names{i}, ' PS']);%scenarios==s
            mean_periodograms{i, indx} = ps;
            hold on;

            %CHANGE THIS BACK TO 'a'
            if strcmp(s, 'a') || strcmp(s, 'cmip6_n')
                unforced_mmms = A.piC_resampled_bootstrapped.r_means(:,:,indx);
                [ps, f, CI] = pmtm(unforced_mmms',[],PAD, 1, 'ConfidenceLevel', .85);
                periods = 1./f(2:end); low = mean(CI(2:end, 1:2:end), 2); high = mean(CI(2:end, 2:2:end), 2);
                %[periods, mp, low, high, ~, ~] = ps_multi(unforced_mmms, obsT, .17,PAD, 2);
                fill([periods;flip(periods)],[low;flip(high)],'y','FaceAlpha', .2 ,'linestyle','-', 'EdgeColor', 'y', 'DisplayName', 'bootstrapped piC MMM PS');
                mean_periodograms{end, indx} = mean(ps, 2);
            end
        end

        %hold on; %p.HandleVisibility = 'off';
    end

    for indx = I
        %subplot(1,I, indx)
        for i=1:length(scenarios)
            s = scenarios{i};
            c = colors{i};
            if(strcmp(s, 'obs'))
                n = [names{i}, ' PS'];
            elseif(subtract_GHG && contains(s, 'g'))
                continue
            else
                n = [names{i}, ' MMM PS'];
            end
            plot(periods, mean_periodograms{i, indx}, '-', 'Color', c, 'DisplayName', n);%scenarios==s
            hold on; 
        end

        legend('location', 'northwest'); xlim([2,end_year-start_year+1]);
        if(strcmp(variable, 'pr'))
            ylabel('Power, ((mm/day)^2/frequency)');
        else
            ylabel('Power, (T^2/frequency)');
        end
        ax = gca; ax.XScale = 'log'; 
        xlabel('Period (Years)'); 
        if(I==1)
            title('Power Spectra of MMMs')
        elseif(subtract_GHG)
            title([A.indices{indx}, ': Forced MMMs - GHG MMM'])
        else
            title([A.indices{indx}, ': Forced MMMs'])
        end
        set(gca, 'XScale', 'linear'); xticks([2, 10:10:114])
        switch realm
            case 'cmip6'
                xlim([2,50]); 
                if(strcmp(variable, 'pr')) 
                    ylim([0,.3]) 
                else
                    ylim([0, 2^(3-indx)]); %hacky, but concise!
                end
            case 'cmip5'
                xlim([2,40]); 
                if(strcmp(variable, 'pr'))
                    ylim([0,.2])
                else
                    ylim([0,.5])
                end
            case 'amip'
                xlim([2,50]);
                ylim([0,3])
        end
    end
end
if(tosave)
    savefig(['figures/', variable, '/', realm, '_Fig3', '.fig']);
    saveas(3, ['figures/', variable, '/', realm, '_Fig3'], 'png');
end

clear mean_periodograms i s c A mmm ps f CI periods fa ls unforced_mmms n ax high low tm indx
%xticks([2,5,10, 20, 50, 100]); ylim([0,.3])
%(automate the tickmarks and switching to a log y-axis)

%TODO MADE it to HERE fixing my code!
%% MAKE FIGURE 5
figure(5); clf;

switch realm
    case 'amip'
        fig5_scenarios = scenarios(1:2);
        fig5_names_1 = {'a. amip-piF simulations', 'b. amip-hist simulations'};
        intense_colors = {[1,.7,0],[0, 127, 0]/255};
        mild_colors = {max(min(intense_colors{1}*2, [1,1,1]), [.4,.4,.6]), max(min(intense_colors{2}*1.5, [1,1,1]), [.4,.6,.4])};
        fnames = names(1:2);
        obs_styles = {'-', '-'};
        %TODO define common_models for amip runs!
    otherwise
        fig5_scenarios = {single_scenario, piC_s};
        fig5_names_1 = {{'Observations and';'ALL simulations'},{'Observations - ALL MMM and';'piC simulations'}};
        mild_colors = {[125,132,241]/255,[255,215,0]/255};
        intense_colors = {'b', [237, 177, 32]/255};
        fnames = {names{1}, 'piC'};
        obs_styles = {'-', '-.'};
end    
AA = load(make_analysis_filename(variable, scenarios{1}, start_year, end_year, '500'));
common_models = AA.indiv.models;

fig5_names_2 = {'c.', 'd.'};
xn = length(fig5_scenarios);
if strcmp(variable, 'ts')
    xn = xn + 1;
    fig5_scenarios = [fig5_scenarios(1), {'noGHG'}, fig5_scenarios(2)];
    fig5_names_1 = [fig5_names_1(1), {{'Observations - GHG MMM and';'ALL Simulations - GHG MMM'}}, fig5_names_1(2)];
    fnames_1 = [fnames(1), {'noGHG'}, fnames(2)];
    obs_styles = [obs_styles, {':'}];
end
if(I==1)
    yn = 2;
else
    yn = I;
end

for basin = 1:I
    yls = nan(1, xn);
    for i = 1:xn
        s = fig5_scenarios{i};
        %this is a little reduntant; but I'm preparing the piC simulations
        %to look like the historical simulations so I can run the scaling
        %there.
        if(strcmp(realm, 'cmip6') && strcmp(variable, 'pr'))
	    h_all = load(make_data_filename('pr', start_month, end_month, 'cmip6_piC','all'));
            h_all.runs(h_all.runs==0)=nan;
            T = table(h_all.model, h_all.runs(:,:,basin), h_all.time, 'VariableNames', {'model', 'runs', 'time'}); 
            T = split_piC_runs(T, end_year - start_year + 1);
            T.time = repmat(start_year:end_year, size(T,1),1);
            h_all = table2struct(T, 'ToScalar', true); clear T;
        elseif ~strcmp(s, 'noGHG')
	    h_all = load(make_data_filename(variable,start_month, end_month,  s, 'all'));
        end
        if(contains(s, 'piC'))
            h_all.runs(h_all.runs==0)=nan;
            T = table(h_all.model, h_all.runs(:,:,basin), h_all.time, 'VariableNames', {'model', 'runs', 'time'}); %pick just NARI
            %use scaling from the historical simulation just before!
            [to_use, loc] = ismember(T.model(:,2), scaling.models); %get rid of models that don't have ALL runs
            T = T(to_use,:);
            T.mMTR = scaling.mMTR(loc(to_use));
            T.color_index_mean = scaling.color_index_mean(loc(to_use));
            clear scaling
            T = split_piC_runs(T, end_year - start_year + 1);
            h_runs = T.runs; h_models = T.model(:,2);
            [scaling.models,idx,scaling.groupings] = unique(h_models, 'stable'); %list of unique models; groupings for each model
            [~,~,scaling.groupings2] = unique(h_models(idx,1)); 
            scaling.mMTR = T.mMTR(idx);
            scaling.color_index_mean = T.color_index_mean(idx);
        elseif(strcmp(s, 'noGHG'))
            %leave scaling as it is. 
            h_runs = h_runs - (G.MMM.MMM(:,:,basin) - mean(G.MMM.MMM(:,:,basin), 2));
        else
            h_runs = h_all.runs(:, ismember(single(h_all.time(1,:,:)), obsT), basin); %pick just NARI
            %+(SCALE-1)*(mmm-mean(mmm, 2));%(:, 1:99)'; %so individual runs will be vertical for PS
            h_models = h_all.model(:,2);
            T = table(h_models, h_all.model(:,1), h_runs, 'VariableNames', {'model', 'institution', 'runs'});
            [to_use, loc] = ismember(T.institution, common_models);
            
            T = T(to_use,:);
            h_runs = T.runs; h_models = T.model;
            %CNRM = contains(h_models, 'CNRM') & strcmp(realm, 'amip');
            scaling = calc_scaling(h_runs, h_models, obs_var(:,:,basin));
        end
        if(strcmp(variable, 'ts'))
            h_all.indices = {'NA', 'GT', 'NARI'};
            scaling.basin_index = h_all.indices{basin};
        end
        scaling.basic_name = fig5_names_1{i};
        FACTOR = mean(obs_var(:,:,basin),2)./scaling.mMTR(scaling.groupings); 
        mean(FACTOR)
        %TODO: remove sqrt
        normalized_runs = (h_runs).*FACTOR;
        normalized_runs_anomalies = normalized_runs - mean(normalized_runs,2);

        if(strcmp(variable, 'ts'))
            subplot(xn, yn, basin+(i-1)*yn)
        else
            subplot(yn,xn,i+(basin-1)*xn)
        end
        if(contains(s, 'piC'))
	     H = load(make_analysis_filename(variable, single_scenario, start_year, end_year, '500'));
            obs_anom_res = obs_var_a(:,:,basin) - (H.MMM.MMM(:,1:length(obs_var_a(:,:,basin)),basin) - mean(H.MMM.MMM(:,1:length(obs_var_a(:,:,basin)),basin),2));
        elseif(strcmp(s, 'noGHG'))
            obs_anom_res = obs_var_a(:,:,basin) - (G.MMM.MMM(:,1:length(obs_var_a(:,:,basin)),basin) - mean(G.MMM.MMM(:,1:length(obs_var_a(:,:,basin)),basin),2));
        else
            obs_anom_res = obs_var_a(:,:,basin);
        end
        %obs_anom_res = obs_anom_res(:,:,basin);%pick just NARI
        %TODO
        [yl, periods, periodograms, ~, ps, CI] = plot_all(h_runs, map, scaling, obs_anom_res, obsT, variable, obs_styles{i}); 
        if(strcmp(variable, 'pr'))
            mean_periodogram = mean(splitapply(@(X) mean(X,1), periodograms, scaling.groupings2),1);
            plot(periods', mean_periodogram(2:end), '--', 'Color', intense_colors{i}, 'DisplayName', ['Mean of ', s, ' PS']);
        end
        set(gca, 'FontSize', 12); 
        yls(i) = yl(2);

        if(I==1)
            subplot(yn,xn,i+xn);
            plot(periods, ps, ['k',obs_styles{i}], 'DisplayName', 'Observed Spectrum', 'Linewidth', 2); hold on;
            fill([periods;flip(periods)],[CI(2:end,1);flip(CI(2:end,2))],'k','FaceAlpha', .05 ,'linestyle',':', 'EdgeColor', 'k', 'EdgeAlpha', .75,'HandleVisibility', 'off');

            %TODO ps_multi doesn't group by model! Each run counts separately. Do I want that?
            [periods, periodograms, CI, CI2] = pmtm_multi(normalized_runs_anomalies, obsT, .17, PAD, scaling.groupings);
            periods = periods(2:end); 
            fill([periods,flip(periods)],[CI2(1,2:end),flip(CI2(2,2:end))],mild_colors{i},'FaceAlpha', .3 ,'linestyle','none', 'DisplayName', ['95% range ',fnames{i},' spectra by model']);
            fill([periods,flip(periods)],[CI(1,2:end),flip(CI(2,2:end))],intense_colors{i},'FaceAlpha', .3 ,'linestyle','none', 'DisplayName', ['66% range ',fnames{i},' spectra by model']);
            plot(periods, mean(periodograms(:,2:end), 1), '--', 'Color', intense_colors{i}, 'DisplayName', ['Mean of ', fnames{i}, ' spectra']);%scenarios==s

            set(gca, 'FontSize', 12);
            
            figure(5+i); clf; plot_all(normalized_runs, map, scaling, obs_anom_res, obsT, variable, obs_styles{i}); 
            xlim([2,end_year-start_year+1]); xlim([2,50]); figure(5)

            if(strcmp(s, single_scenario))%~any([contains(s, 'piC'), contains(s, 'amip')]))
                for f=2:length(scenarios)
                    s = scenarios{f};
                    if(contains(s, 'fast') || contains(s, 'obs') || strcmp(realm, 'amip'))
                        continue
                    end
                    c = colors{f};%scenarios==s);
		    A = load(make_data_filename(variable, start_month, end_month, s, 'all'));
                    [~,l_idx,local_groupings] = unique(A.model(:,2), 'stable');
                    [~,~,local_groupings2] = unique(A.model(l_idx,1), 'stable');
                    forced_runs = A.runs(:,:,basin); %pick just NARI
                    normalized_forced_runs = forced_runs./mean(forced_runs,2)*mean(obs_var(:,:,basin),2); %I'm not grouping by model for this mean.
                    normalized_forced_runs_anomalies = normalized_forced_runs - mean(normalized_forced_runs,2);
                    [periods, periodograms, ~, ~] = pmtm_multi(normalized_forced_runs_anomalies, obsT, [], PAD, local_groupings, local_groupings2);
                    plot(periods, mean(periodograms, 1), '--', 'Color', c, 'DisplayName', ['Tiered Mean of ', names{f}, ' spectra']);%scenarios==s
                end
            end
            legend('location', 'northwest')
            ylabel('Power, Scaled', 'FontWeight', 'bold'); ylim(yl); 
            xlabel('Period (years)'); xlim([2,50]);
            title(fig5_names_2{i})
        end
        clear A c CI CI2 forced_runs h_CI h_mean_periodogram half i idx
        clear I_idx local_groupings local_groupings2 normalized* periodograms periods
        clear ps s
    end
    for i=1:xn
        if(strcmp(variable, 'ts'))
            subplot(xn, yn, basin+(i-1)*yn)
        else
            subplot(yn,xn,i+(basin-1)*xn)
        end
        ylim([0, min(yls)])
    end
end

if(tosave)
    savefig(['figures/', variable, '/', realm, '_Fig5.fig']); %TODO!
    saveas(5, ['figures/', variable, '/', realm, '_Fig5'], 'png');
end

%{
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

%% Functions
function [scaling] = calc_scaling(h_runs, h_models, obs_var)
    total_rainfall = mean(h_runs, 2); %why did I switch these? TODO
    %scaling.indices = (1:length(total_rainfall))'; %just a number for each element for use in indexing later
    [scaling.models,idx,scaling.groupings] = unique(h_models, 'stable'); %list of unique models; groupings for each model
    [~,~,scaling.groupings2] = unique(h_models(idx,1)); 
    scaling.mMTR = splitapply(@mean, total_rainfall, scaling.groupings);
    half = max(abs(scaling.mMTR-mean(obs_var,2)));
    scaling.color_index_mean = floor(1+254/2*(1+(scaling.mMTR-mean(obs_var,2))/half));
end

function [yl, periods, periodograms, h_CI, obs_ps, CI] = plot_all(m_runs, color_map, scaling, obsv, obsT, variable, obs_style)
    if nargin < 6
        obs_style = '-';
    end
    global PAD
    %I should smooth these as well to make up for the smaller means
    [periods, periodograms, ~, h_CI] = pmtm_multi(m_runs, obsT, [], PAD, scaling.groupings);
    for i = unique(scaling.groupings)'
        plot(periods, periodograms(i,:), 'color', color_map(scaling.color_index_mean(i),:), 'Linewidth', 1, 'HandleVisibility', 'off','DisplayName', scaling.models{i})
        hold on;
    end
    [obs_ps, f, CI] = pmtm(obsv',[],PAD, 1, 'ConfidenceLevel', .95);
    obs_ps = obs_ps(2:end); periods = 1./f(2:end);
    plot(periods, obs_ps, obs_style, 'color', 'k', 'DisplayName', 'Observations', 'Linewidth', 2)
    fill([periods;flip(periods)],[CI(2:end,1);flip(CI(2:end,2))],'k','FaceAlpha', .05 ,'linestyle',':', 'EdgeColor', 'k', 'EdgeAlpha', .75,'HandleVisibility', 'off');
    if(strcmp(variable, 'ts'))
        ylabel(scaling.basic_name, 'FontWeight', 'bold'); 
        title(scaling.basin_index, 'FontWeight', 'bold');
    else
        ylabel('Power, True Values', 'FontWeight', 'bold'); 
        title(scaling.basic_name, 'FontWeight', 'bold'); 
    end
    xlabel(''); xlim([2,50]); yl = ylim;
end

function [T] = split_piC_runs(T, LENGTH)
    reps = floor(sum(~isnan(T.runs),2)/LENGTH);
    T = repelem(T, reps, 1); reps = reps(reps~=0);
    count = @(x) {(1:x)'};
    groups = splitapply(count, reps, (1:length(reps))');
    groups = cat(1,groups{:});
    runs = splitapply(@(x1,x2){circshift(x1,x2)}, table(T.runs, -(groups-1)*LENGTH), (1:length(groups))');
    time = splitapply(@(x1,x2){circshift(x1,x2)}, table(T.time, -(groups-1)*LENGTH), (1:length(groups))');
    T.runs = cat(1,runs{:}); T.time = cat(1,time{:});
    T.runs = T.runs(:,1:LENGTH); T.time = T.time(:,1:LENGTH);
end

%NA: CNRM-ESM2-1 p1 (pink), IPSL-CM6A-LR p1 (blue), CNRM-CM6-1 p1 (grey)
%NARI: CNRM-ESM2-1 p1 (pink), CNRM-CM6-1 p1 (grey pink)

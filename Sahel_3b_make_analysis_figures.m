mic=false;

scenarios = {'r'};%'v'};%'a6'};%'e'};%'h','a','n','g'};
single_scenario='r';%'v';%'a6';%'e';%'amip';%'h';%
short = false;%true; %TODO I've totally bastardized what this was supposed to mean -- redo later
single_name = 'AMIP+RAD';%'Vanilla AMIP';%'CMIP6-AMIP';%'ERA-20CM';%'ALL';%'AMIP';%
month0 = 'Jul';%, "Jun", "Sep"];
month1 = 'Sep';%, "Jul", "Oct"];
start_year = 1901;

dts = [""];%"detrended",    ,            This is implemented
fls = ["last"];%, "first"       ,        This is implemented
colors = 'bmrg';%[blue; magenta; red; green];
names = [{'AMIP+RAD'}];%'Vanilla AMIP'}];%[{'CMIP6-AMIP'}]; %[{'ERA-20CM'}]; %[{'ALL'}, {'AA'}, {'NAT'}, {'GHG'}];

N = 500;

load('model_output/historical_precipitation.mat')

%% MAKE SECOND HALF OF FIGURE 1
close all;
openfig([single_scenario, '_Fig1.fig']);
hold on; %hold off; clf; 
aname = ['analysis/',single_scenario,'_N', num2str(N)]; %TODO
if(short)
    aname = [aname, '_', num2str(start_year)];
end
A = load([aname, '.mat']);
r = A.MMM.r; r_indiv = A.indiv.r; r_bootstrap = A.historical_bootstrapped.rs; r_r = A.indiv_runs.r; 
e = A.MMM.e; e_indiv = A.indiv.e; e_bootstrap = A.historical_bootstrapped.es; e_r = A.indiv_runs.e;
%indiv_plot(r_r, 3); hold on; indiv_plot(e_r, 4); hold on;
indiv_plot(r_indiv, 3); hold on; analysis_plot(r, r_bootstrap, 'b', 'MMM', 1, 2, 3);
indiv_plot(e_indiv, 4); hold on; analysis_plot(e, e_bootstrap, 'b', 'MMM', 1, 2, 4);
finishfig(1,2,3,'c. Correlation with 20C Observations', '', 0, 'northwest'); 
finishfig(1,2,4,'d. RMSE with 20C Observations', 'Fraction of Observed Variance', 1, 'northeast');
savefig([single_scenario, '_Fig1_all.fig']);
%{
% % MAKE FIGURE 4
if(~mic)
    F = 4; figure(F); hold off; clf;
    for j = 1:length(scenarios)
        scenario = scenarios{j};
        color = colors(j);
        name = names{j};
        for k = 1:length(fls)
            fl = char(fls(k));
            for l = 1:length(dts)
                dt = char(dts(l));

                aname = ['analysis/', scenario, '_N', num2str(N)];
                if(short)
                    aname = [aname, '_', num2str(start_year)];
	        end
		A = load([aname, '.mat']);

                r = A.MMM.r; r_indiv = A.indiv.r; r_bootstrap = A.historical_bootstrapped.rs; 
                e = A.MMM.e; e_indiv = A.indiv.e; e_bootstrap = A.historical_bootstrapped.es; 

                analysis_plot(r, r_bootstrap, color, name, F, 1, 1);
                analysis_plot(e, e_bootstrap, color, name, F, 1, 2); 
%    {
                r_resample = A.piC_resampled_bootstrapped.rs; e_resample = A.piC_resampled_bootstrapped.es; 
                analysis_plot(r, r_resample, color, 'Pre-Industrial Control', F, 1, 1)
                analysis_plot(e, e_resample, color, 'Pre-Industrial Control', F, 1, 2)
%  }
            end
        end
    end
    finishfig(F,1,1,'a. Correlation with 20C Observations', '', 0); 
    finishfig(F,1,2,'b. Scaled RMSE with 20C Observations', '', 1);
    savefig([scenario, '_Fig4.fig']);
end
%}
%% MAKE MIC's figures
if(mic)
    F = 20; figure(F); clf; hold off; clf;

    aname = ['analysis/', single_scenario, '_N', num2str(N)];%
    if(strcmp(single_scenario, 'amip'))
        aname=[aname,'.mat'];
    elseif(short)
        aname=[aname, '_', num2str(start_year), '.mat'];% '_short.mat'];
    else
        aname=[aname,'.mat'];
    end
    A = load(aname);

    std_O = std(prcp); 
    [rs_nf, as_uncorr, Bs] = calc_Mic(std_O, A.historical_bootstrapped.b_means, A.historical_bootstrapped.rs);
    [r_nf, a_uncorr, B] = calc_Mic(std_O, A.MMM.MMM, A.MMM.r);
    %{
    [~, a_piC, ~] = calc_Mic(std_O, A.MMM.MMM, A.piC_resampled_bootstrapped.rs);

    piC = load('model_output/piC_all.mat'); piC.runs(piC.runs==0)=NaN;
    piC_runs = piC_select(N,length(A.MMM.MMM), piC.runs);
    piC_rs = corr(A.MMM.MMM', piC_runs');
    piC_Bs_sq = 1-var(piC_runs,0,2)/var(prcp);
    piC_B_sq = 1-mean(var(piC_runs,0,2))/var(prcp);

    piC_Bs_a_sq = 1-a_uncorr^2.*var(piC_runs,0,2)/var(prcp);
    piC_B_a_sq = 1-a_uncorr^2*mean(var(piC_runs,0,2))/var(prcp);
    %}

    analysis_plot(r_nf, rs_nf, 'b', [single_name,' MMM'], F, 1, 1);
    %analysis_plot(r_nf, piC_rs', 'y', 'piC runs', F, 1, 1);
    analysis_plot(a_uncorr, as_uncorr, 'b', [single_name,' MMM'], F, 1, 2); 
    %analysis_plot(a_uncorr, a_piC, 'y', 'piC MMMs', F, 1, 2); 

    finishfig(F,1,1,'a. Correlation between Noise and Forced Signal if \alpha = 1', '', 0); 
    finishfig(F,1,2,'b. \alpha (Scaling of Forced Signal) for Uncorrelated Noise', '', 1);
    legend('location', 'northwest');
    savefig([single_scenario, '_', num2str(start_year), '_MicFig.fig']);

    figure(F+1); clf; hold on;
    analysis_plot(B^2, Bs.^2, 'b', 'ALL MMM', F+1, 1, 1);
    %analysis_plot(piC_B_sq, piC_Bs_sq, 'y', 'piC runs', F+1, 1, 1);
    %analysis_plot(piC_B_a_sq, piC_Bs_a_sq, 'g', 'alpha * piC runs', F+1, 1, 1);
    finishfig(F+1,1,1,'\beta^2 (Ratio of the Forced to Observed Variance) according to the...', '', 0); 
    legend('location', 'northeast');
end

%% Functions
function [runs] = piC_select(N, T, GMs)
    s = size(GMs); num_models=s(1);
    N_times = floor(N/num_models); N=N_times*num_models;
    idx_raw = cumsum(ones(s), 2);
    length_run = sum(T);
    length_of_GMs = sum(~isnan(GMs),2);  
    runs = nan(N, length_run);
    for i=1:N_times
        rand_end_offset = nan(num_models, 1);
        for mn = 1:num_models
            rand_end_offset(mn) = randi([length_run, length_of_GMs(mn)]);
        end 
        idx = (idx_raw > rand_end_offset - length_run) & (idx_raw <= rand_end_offset);
        m_t = GMs'; idx_t = idx';
        rng = (1:num_models)+(i-1)*num_models;
        runs_t = vec2mat(m_t(idx_t), length_run);
        runs_t = runs_t - mean(runs_t, 2);
        runs(rng, :) = runs_t;
    end
end

function[r_nf, a, Bs] = calc_Mic(std_O, means, rs)
    stds_F = std(means, 0, 2); Bs = stds_F/std_O; 
    a = rs./Bs; r_nf = (rs-Bs)./(1+Bs.^2-2*rs.*Bs);
end

function[] = indiv_plot(indiv, sp_num)
    colorList = get(gca,'ColorOrder');
    cyan = colorList(6,:);
    subplot(2,2,sp_num);
    histogram(indiv, 'Normalization', 'pdf',...
        'FaceColor', cyan,...
        'DisplayName', 'IMs',...
        ...%'BinWidth', .05);
        'NumBins', 7);
    hold on;
end

function[] = analysis_plot(v, bootstrap, color, name, fig_num, sp_x, sp_num)%indiv)
    figure(10)
    h = histfit(bootstrap, [], 'kernel');
    %TODO I won't always want this number of subplots...
    figure(fig_num);
    if(fig_num ~= 21); subplot(sp_x,2,sp_num); end %oyoyoy
    lineplot = h(2); x = lineplot.XData; y = lineplot.YData; y = y / ((x(2) - x(1)) * sum(y));
    idx_v = find(abs(x - v) == min(abs(x - v)));
    if(contains(name, 'piC')||contains(name, 'Control'))
        plot(x, y, [color, ':'], 'LineWidth', 2, 'DisplayName', name); hold on;
        [low, high] = confidence_interval(bootstrap, .05); 
        if(v > median(bootstrap))
            plot(high*[1,1], get(gca, 'ylim'), [color, '--'], 'HandleVisibility', 'off');
        else
            plot(low*[1,1], get(gca, 'ylim'), [color, '--'], 'HandleVisibility', 'off');
        end
    else
        plot(x, y, [color, '-*'], 'MarkerIndices', [idx_v], 'LineWidth', 2, 'DisplayName', [name, ' = ', num2str(v)]);
    end
    hold on;
end

function[] = finishfig(fig_num, sp_x, sp_num, t, xl, s, l)
    if(nargin < 7)
        l = 0;
    end
    figure(fig_num); 
    if(fig_num ~= 21); subplot(sp_x,2,sp_num); end %TODO 
    set(gca,'FontSize',16); ylabel('Density'); 
    yl = ylim;
    plot(s*[1,1], yl, 'k-', 'HandleVisibility','off');
    if(nargin >= 7), location = l; elseif(s), location = 'northeast'; else, location = 'northwest'; end
    legend('show', 'Location', location)
    title(t); ylim(yl); xlabel(xl);
end

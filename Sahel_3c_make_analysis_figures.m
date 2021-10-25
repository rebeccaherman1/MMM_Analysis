%For TS just looks at NARI.
mic=false;
tosave = false;

realm = 'cmip6_fast';
variable = 'pr';
start_year = 1901; end_year = 2014;
short = false;
start_month = 7;
end_month = 9;

switch realm
    case 'cmip6'
        scenarios = {'cmip6_h', 'cmip6_a', 'cmip6_n', 'cmip6_g', 'amip-hist', 'amip-piF'};% The Fast component is SUPER WIDE! 'cmip6_fast'};
        colors = {'b', 'm', [0.60,0.20,0.00], 'g', [0, 127, 0]/255, [1,.7,0],[126, 47, 142]/255};
        styles = {'-', '-', '-', '-', '-', '-'};
        names = {'ALL 6', 'AA 6', 'NAT 6', 'GHG 6', 'amip-hist', 'amip-piF', 'Implied Fast Component'};
        mdgnd = [0,1,1]; %cyan
        %I COULD automate adding cmip5, but I don't feel like it...
        if strcmp(variable, 'ts')
            scenarios = scenarios(1:4); %the others will get automatically shortened. 
        end
        if(short)
            scenarios = [scenarios(1:4), {'h', 'a', 'n', 'g'}];
            colors = [colors(1:4), {[0.00,0.45,0.74], [0.75,0.00,0.75], 'r', [0.47,0.67,0.19]}];
            styles = [styles(1:4), {'-.','-.','-.','-.'}];
            names = [names(1:4), {'ALL', 'AA', 'VA', 'GHG'}];
            end_year = 2003;
        end
    case 'cmip6_fast'
        scenarios = {'cmip6_hfast', 'cmip6_afast', 'cmip6_nfast', 'cmip6_gfast'};% The Fast component is SUPER WIDE! 'cmip6_fast'};
        colors = {'b', 'm', [0.60,0.20,0.00], 'g', [0, 127, 0]/255, [1,.7,0],[126, 47, 142]/255};
        styles = {'-', '-', '-', '-', '-', '-'};
        names = {'ALL fast', 'AA fast', 'NAT fast', 'GHG fast', 'amip-hist', 'amip-piF', 'Implied Fast Component'};
        mdgnd = [0,1,1]; %cyan
    case 'cmip5'
        scenarios = {'h', 'a', 'n', 'g'};
        colors = {'b', 'm', [0.60,0.20,0.00], 'g'};
        styles = {'-', '-', '-', '-'};
        names = {'ALL', 'AA', 'VA', 'GHG'};
        mdgnd = [0,1,1]; %cyan
        end_year = 2003;
    case 'r'
        %scenarios = {'cmip6_r', 'v'};
        scenarios = {'amip-hist'};
        styles = {'-'};
        colors = {[0, 127, 0]/255};
        names = {'amip-hist'};
        mdgnd = max(min(colors{1}*2, [1,1,1]), [.4,.4,.6]);
    case 'v'
        scenarios = {'amip-piF'};
        colors = {[1,.7,0]};
        styles = {'-'};
        names = {'amip-piF'};
        mdgnd = max(min(colors{1}*1.5, [1,1,1]), [.4,.6,.4]);
    otherwise
        %a6? e? p? amip? (that last one's cmip5)
end
single_scenario = scenarios{1};

dts = [""];%"detrended",    ,            This is implemented
fls = ["last"];%, "first"       ,        This is implemented

N = 500;

%load('data/', variable, '/observations.mat')

%% MAKE SECOND HALF OF FIGURE 1

%TODO add NARI for amip figures.
%TODO do I want to see fig 1 for SST? or the other figures are enough?

close all;
 %hold off; clf; 

if(strcmp(variable, 'pr') && ~short)
    aname = make_analysis_filename(variable, start_month, end_month, single_scenario,start_year, end_year, N);
    A = load(aname);
    openfig(['figures/', variable, '/', single_scenario, '_Fig1p_', num2str(start_year), '-', num2str(end_year), '.fig']);
    hold on;

    r = A.MMM.r; r_indiv = A.indiv.r; r_bootstrap = A.historical_bootstrapped.rs; r_r = A.indiv_runs.r; 
    e = A.MMM.e; e_indiv = A.indiv.e; e_bootstrap = A.historical_bootstrapped.es; e_r = A.indiv_runs.e;
    %indiv_plot(r_r, 3); hold on; indiv_plot(e_r, 4); hold on;
    indiv_plot(r_indiv, 3, mdgnd); hold on; analysis_plot(r, r_bootstrap, colors{1}, styles{1}, 'MMM', 1, 2, 3);
    indiv_plot(e_indiv, 4, mdgnd); hold on; analysis_plot(e, e_bootstrap, colors{1}, styles{1}, 'MMM', 1, 2, 4);
    finishfig(1,2,3,'c. Correlation with 20C Observations', '', 0, 'northwest'); 
    finishfig(1,2,4,'d. RMSE with 20C Observations', 'Fraction of Observed Variance', 1, 'northeast');
    if(tosave)
        savefig(['figures/', variable, '/', single_scenario, '_Fig1_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.fig']);
        saveas(1, ['figures/', variable, '/', single_scenario, '_Fig1_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N)], 'png');
    end
end

%% MAKE FIGURE 4
%TODO if I want to be able to make these figures for TS, I have to make
%multiple subplots and add a for loop over I. 

%Just NARI.
if(contains(realm, 'cmip'))
    F = 4; figure(F); hold off; clf;
    for j = 1:length(scenarios)
        scenario = scenarios{j};
        color = colors{j};
        style = styles{j};
        name = names{j};
        for k = 1:length(fls)
            fl = char(fls(k));
            for l = 1:length(dts)
                dt = char(dts(l));

		aname = make_analysis_filename(variable, scenario, start_year, end_year, N);
                A = load(aname);

                r = A.MMM.r; r_bootstrap = A.historical_bootstrapped.rs; 
                e = A.MMM.e; e_bootstrap = A.historical_bootstrapped.es; 
                
                if(strcmp(variable, 'ts'))
                    r = r(:,:,3); r_bootstrap = r_bootstrap(:,:,3);
                    e = e(:,:,3); e_bootstrap = e_bootstrap(:,:,3);
                end

                analysis_plot(r, r_bootstrap, color, style, name, F, 1, 1);
                analysis_plot(e, e_bootstrap, color, style, name, F, 1, 2); 
                
                if(isfield(A, 'piC_resampled_bootstrapped')) %UMMMM why does amip-piF have piC??? it shouldn't....
                    r_resample = A.piC_resampled_bootstrapped.rs; e_resample = A.piC_resampled_bootstrapped.es; 
                    if(strcmp(variable, 'ts'))
                        r_resample = r_resample(:,:,3); e_resample = e_resample(:,:,3);
                    end
                    analysis_plot(r, r_resample, color, style, 'Pre-Industrial Control', F, 1, 1)
                    analysis_plot(e, e_resample, color, style, 'Pre-Industrial Control', F, 1, 2)
                end
            end
        end
    end
    finishfig(F,1,1,'a. Correlation with 20C Observations', '', 0); 
    finishfig(F,1,2,'b. sRMSE with 20C Observations', '', 1);
    if(tosave)
        savefig(['figures/', variable, '/', realm, '_Fig4', '_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.fig']);
    end
end

%% MAKE MIC's figures
%{
if(true)%mic)
    F = 20; figure(F); clf; hold off; clf;

    if(short && ~strcmp(single_scenario, 'amip'))
        aname = ['analysis/',single_scenario,'_N', num2str(N), '_', num2str(start_year)];
    elseif(extended)
        aname = ['analysis/',single_scenario,'_N', num2str(N), '_extended']; %TODO
    else
        aname = ['analysis/',single_scenario,'_N', num2str(N)]; %TODO
    end
    A = load([aname, '.mat']);

    std_O = std(prcp); 
    [rs_nf, as_uncorr, Bs] = calc_Mic(std_O, A.historical_bootstrapped.b_means, A.historical_bootstrapped.rs);
    [r_nf, a_uncorr, B] = calc_Mic(std_O, A.MMM.MMM, A.MMM.r);
    
    [~, a_piC, ~] = calc_Mic(std_O, A.MMM.MMM, A.piC_resampled_bootstrapped.rs);

    piC = load('data/cmip6_piC_all.mat'); piC.runs(piC.runs==0)=NaN;
    piC_runs = piC_select(N,length(A.MMM.MMM), piC.runs);
    piC_rs = corr(A.MMM.MMM', piC_runs');
    piC_Bs_sq = 1-var(piC_runs,0,2)/var(prcp);
    piC_B_sq = 1-mean(var(piC_runs,0,2))/var(prcp);

    piC_Bs_a_sq = 1-a_uncorr^2.*var(piC_runs,0,2)/var(prcp);
    piC_B_a_sq = 1-a_uncorr^2*mean(var(piC_runs,0,2))/var(prcp);
    

    analysis_plot(r_nf, rs_nf, 'b', [single_name,' MMM'], F, 1, 1);
    analysis_plot(r_nf, piC_rs', 'y', 'piC runs', F, 1, 1);
    analysis_plot(a_uncorr, as_uncorr, 'b', [single_name,' MMM'], F, 1, 2); 
    %analysis_plot(a_uncorr, a_piC, 'y', 'piC MMMs', F, 1, 2); 

    finishfig(F,1,1,'a. Correlation between Noise and Forced Signal if \alpha = 1', '', 0); 
    finishfig(F,1,2,'b. \alpha (Scaling of Forced Signal) for Uncorrelated Noise', '', 1);
    legend('location', 'northwest');
    if(tosave)
        savefig([single_scenario, '_MicFig.fig']);
    end
    
    figure(F+1); clf; hold on;
    analysis_plot(B^2, Bs.^2, 'b', 'ALL MMM', F+1, 1, 1);
    %analysis_plot(piC_B_sq, piC_Bs_sq, 'y', 'piC runs', F+1, 1, 1);
    %analysis_plot(piC_B_a_sq, piC_Bs_a_sq, 'g', 'alpha * piC runs', F+1, 1, 1);
    finishfig(F+1,1,1,'\beta^2 (Ratio of the Forced to Observed Variance) according to the...', '', 0); 
    legend('location', 'northeast');
end
%}
%% Functions
function [runs] = piC_select(N, T, GMs)
    length_run = sum(T);    
    length_of_GMs = sum(~isnan(GMs),2);
    too_short = length_of_GMs<length_run;
    
    length_of_GMs = length_of_GMs(~too_short);
    GMs = GMs(~too_short,:); 

    s = size(GMs); num_models=s(1);
    N_times = floor(N/num_models); N=N_times*num_models;
    idx_raw = cumsum(ones(s), 2);
    
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

function[] = indiv_plot(indiv, sp_num, clr)
    subplot(2,2,sp_num);
    histogram(indiv, 'Normalization', 'pdf',...
        'FaceColor', clr,...
        'DisplayName', 'IMs',...
        ...%'BinWidth', .05);
        'NumBins', 7);
    hold on;
end

function[] = analysis_plot(v, bootstrap, color, style, name, fig_num, sp_x, sp_num)%indiv)
    figure(10)
    h = histfit(bootstrap, [], 'kernel');
    %TODO I won't always want this number of subplots...
    figure(fig_num);
    if(fig_num ~= 21); subplot(sp_x,2,sp_num); end %oyoyoy When TH do I put 21?
    lineplot = h(2); x = lineplot.XData; y = lineplot.YData; y = y / ((x(2) - x(1)) * sum(y));
    idx_v = find(abs(x - v) == min(abs(x - v)));
    if(contains(name, 'piC')||contains(name, 'Control'))
        plot(x, y, ':', 'Color', color, 'LineWidth', 2, 'DisplayName', name); hold on;
        [low, high] = confidence_interval(bootstrap, 1, .001); 
        if(v > median(bootstrap))
            plot(high*[1,1], get(gca, 'ylim'), '--', 'Color', color, 'HandleVisibility', 'off');
        else
            plot(low*[1,1], get(gca, 'ylim'), '--', 'Color', color, 'HandleVisibility', 'off');
        end
    else
        plot(x, y, [style, '*'], 'Color', color, 'MarkerIndices', idx_v, 'LineWidth', 2, 'DisplayName', [name, ' = ', num2str(v, '%5.2f')]);
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

historical = true;
pC = true;
save = false;
standardize = false;

scenarios = {'cmip6_h','cmip6_a','cmip6_n','cmip6_g'};%'amip'};%,; 
scenario_names = {'ALL 6', 'AA 6', 'NAT 6', 'GHG 6'};
month0 = 'Jul';%, "Jun", "Sep"];
month1 = 'Sep';%, "Jul", "Oct"];

%dts = ["detrended", ""];%   ,            This is implemented
%fls = ["last", "first"];%       ,        This is implemented
%dis = [""];% "Decadal", "Interannual"  TODO!
fl = '';%'last';
N=500;

scenario_colors = 'bmrgc';

global ref_T_years; 
start_year = 1901;%1979; 
end_year=2003;

obs = load('data/historical_precipitation.mat');
timeframe = (obs.T >= start_year & obs.T <= end_year);
ref_T_years = obs.T(timeframe);
prcp_unstandardized = obs.prcp(timeframe);
prcp_anomaly = prcp_unstandardized - mean(prcp_unstandardized);
prcp_std = std(prcp_unstandardized);
prcp_standardized = prcp_anomaly/prcp_std;

figure(2); hold off; clf; hold on; 
%set(gca,'FontSize',12);%15);
%set(gca,'LineWidth',2);
if(standardize)
    fig_title = ['Standardized Series Comparison: ', month0, '-', month1];
    ylabel("Precipitation Anomaly");
else
    fig_title = ['Series Comparison: ', month0, '-', month1];
    ylabel("Precipitation Anomaly (kg/m^2/s)");
end
title(fig_title);

N_indiv = nan(1,length(scenarios));
for j = 1:length(scenarios)
    scenario = char(scenarios(j));
    fprintf("Accessing historical scenario %s\n", scenario);
    GM = load(['data/', scenario, '_GM.mat']);
    N_indiv(j) = length(GM.models);
    clear GM;

    A = load(['Analysis/', scenario, '_N', num2str(N), '.mat']);
    %which PC do we want?
    subplot(2, 2, j); 
    if (standardize)
        p_actual = plot(ref_T_years, prcp_standardized, 'k');
    else
        p_actual = plot(ref_T_years, prcp_anomaly, 'k-');%, 'LineWidth', 2);
    end
    hold on; set(gca,'FontSize',15); 
    plot(ref_T_years, zeros(size(ref_T_years)), 'k--'); %haha I don't think anyone actually ever sees this...
    yyaxis right
    set(gca, 'ycolor', scenario_colors(j))
    ylim([-.375,.375]) %TODO comment out for AMIP
    %pC
    if(pC)
        if(strcmp(fl, 'last'))
            down = A.piC_last.low;
            up = A.piC_last.high;
            scale = mean(std(A.piC_last.r_means, 0, 2));         
        else
            down = A.piC_resampled_bootstrapped.low;
            up   = A.piC_resampled_bootstrapped.high;
            scale = mean(std(A.piC_resampled_bootstrapped.r_means, 0, 2));
        end
        if(standardize)
            down = down/scale;
            up = up/scale;
            ylim([-5,5])
        end

        p_stds = fill([ref_T_years';flipud(ref_T_years')],[down';flipud(up')],'y','FaceAlpha', .3 ,'linestyle','none');
        p_stds.HandleVisibility = 'off';
        %title([scenario, ' piC, r = ', num2str(r), ', rmse = ', num2str(rmsd), ', N = ', num2str(Analysis.n_pC_means)]);
    end

    if(historical)
        mmm = A.MMM.MMM; mmm = mmm-mean(mmm); r = A.MMM.r; rmsd = A.MMM.e;
        up   = A.historical_bootstrapped.high;
        down = A.historical_bootstrapped.low;
        if(standardize)
            scale = std(mmm);
            mmm = mmm/scale;
            down = down/scale;
            up = up/scale;
            ylim([-5,5])
        end

        hold on;
        p_stds = fill([ref_T_years';flipud(ref_T_years')],[down';flipud(up')],scenario_colors(j),'FaceAlpha', .3 ,'linestyle','none');
        p_stds.HandleVisibility = 'off';
        p_mmm = plot(ref_T_years,mmm,[scenario_colors(j),'-']);%, 'LineWidth', 2);%, 'DisplayName', );
        title([scenario, ', r = ', num2str(r), ', rmse = ', num2str(rmsd), ', N = ', num2str(A.historical_bootstrapped.num_models)]);
        xlim([start_year, end_year]);
    end
    title([scenario_names{j}, ', N=', num2str(N_indiv(j))], 'color', scenario_colors(j))
end
subplot(2,2,1)
yyaxis left
ylabel('Precipitation Anomaly (mm/day)')
subplot(2,2,3)
yyaxis left
ylabel('Precipitation Anomaly (mm/day)')
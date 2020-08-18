historical = true;
pC = true;
tosave = false;
start_year = 1901;
end_year = 2014;
variable = 'ts';
realm = 'cmip6';

scenario_colors = 'bmrgc';
long_colors = {'blue', 'magenta', 'red', 'green', 'cyan'};
if(strcmp(realm, 'cmip6'))
    scenarios = {'cmip6_h','cmip6_a','cmip6_n','cmip6_g'};%'amip'};%,; 
    scenario_names = {'ALL 6', 'AA 6', 'NAT 6', 'GHG 6'};
elseif(strcmp(realm, 'cmip5'))
    scenarios = {'h', 'a', 'n', 'g'};
    scenario_names = {'ALL', 'AA', 'NAT', 'GHG'};
else
    scenarios = {realm};
end
zoom = strcmp(variable, 'pr');

fl = '';%'last';
N=500;

global ref_T_years; 

obs = load(['data/', variable, '/observations.mat']);

timeframe = (obs.T >= start_year & obs.T <= end_year);
ref_T_years = obs.T(timeframe)';
start_year = max(start_year, ref_T_years(1));
end_year = min(end_year, ref_T_years(end));
var_unstandardized = obs.var(:,timeframe,:);
var_anomaly = var_unstandardized - mean(var_unstandardized,2);
var_std = std(var_unstandardized,0,2);
var_standardized = var_anomaly./var_std;

figure(2); hold off; clf; hold on; 
g_GM = load(['data/', variable, '/', scenarios{contains(scenarios, 'g')}, '_GM.mat']);
no_ghg = var_anomaly - (g_GM.MMM - mean(g_GM.MMM, 2));

N_indiv = nan(1,length(scenarios));
for j = 1:length(scenarios)
    scenario = char(scenarios(j));
    fprintf("Accessing historical scenario %s\n", scenario);
    GM = load(['data/', variable, '/', scenario, '_GM.mat']);
    A = load(['Analysis/', variable, '/', scenario, '_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.mat']);
    N_indiv(j) = length(GM.models);
    clear GM;
    [~,~,I] = size(A.MMM.MMM);

    if(pC)
        if(strcmp(fl, 'last'))
            pC_down = A.piC_last.low;
            pC_up = A.piC_last.high;
            pC_scale = mean(std(A.piC_last.r_means, 0, 2));         
        else
            pC_down = A.piC_resampled_bootstrapped.low;
            pC_up   = A.piC_resampled_bootstrapped.high;
            pC_scale = mean(std(A.piC_resampled_bootstrapped.r_means, 0, 2));
        end
    end
    if(historical)
        mmm = A.MMM.MMM; mmm = mmm-mean(mmm,2); r = A.MMM.r; rmsd = A.MMM.e;
        up   = A.historical_bootstrapped.high;
        down = A.historical_bootstrapped.low;
    end
    
    %which PC do we want?
    for i=1:I
        if(~contains(realm, 'cmip'))
            yn = 1; xn = 1;
        elseif(I==1)
            %only works if perfect square!
            yn = length(scenarios)^.5; xn = yn;
        	subplot(yn, xn, j); 
        else
            yn = 4; xn = I;
            subplot(yn,xn, I*(j-1)+i)
        end
        plot(ref_T_years, zeros(size(ref_T_years)), 'k--'); 
        hold on; set(gca,'FontSize',15); 
        if(I==1)
            title([scenario_names{j}, ', N=', num2str(N_indiv(j))], 'color', scenario_colors(j))
        else
            ttl = '';
            if(i==1)
                ttl = ['\color{', long_colors{j}, '}', scenario_names{j}, ', N=', num2str(N_indiv(j))];
            end
            if(j==1)
                ttl = [ttl, ' \color{black}', A.indices{i}];
            end
            title(ttl)
        end 
        if(strcmp(variable, 'ts') && (contains(scenario, 'a') || contains(scenario, 'n')) && i<3)
            p_actual = plot(ref_T_years, no_ghg(:,:,i), '-', 'color', [.5 0 .6]);%, 'LineWidth', 2);
        else
            p_actual = plot(ref_T_years, var_anomaly(:,:,i), 'k-');%, 'LineWidth', 2);
        end
        if(zoom)
            yyaxis right
            set(gca, 'ycolor', scenario_colors(j))
            ylim([-.375,.375])
        else
            ylim([-1,1])
        end
        %pC
        if(pC)
            p_stds = fill([ref_T_years';flipud(ref_T_years')],[pC_down(:,:,i)';flipud(pC_up(:,:,i)')],'y','FaceAlpha', .3 ,'linestyle','none');
            p_stds.HandleVisibility = 'off';
            %title([scenario, ' piC, r = ', num2str(r), ', rmse = ', num2str(rmsd), ', N = ', num2str(Analysis.n_pC_means)]);
        end
        if(historical)
            hold on;
            p_stds = fill([ref_T_years';flipud(ref_T_years')],[down(:,:,i)';flipud(up(:,:,i)')],scenario_colors(j),'FaceAlpha', .3 ,'linestyle','none');
            p_stds.HandleVisibility = 'off';
            p_mmm = plot(ref_T_years,mmm(:,:,i),[scenario_colors(j),'-']);%, 'LineWidth', 2);%, 'DisplayName', );
            %title([scenario, ', r = ', num2str(r), ', rmse = ', num2str(rmsd), ', N = ', num2str(A.historical_bootstrapped.num_models)]);
            xlim([start_year, end_year]);
        end
        %still need this for pr... and it has to be at the end, it would
        %seem.
    end
end
if(strcmp(variable, 'pr'))
    lbl = 'Precipitation Anomaly (mm/day)';
else
    lbl = 'SST Anomaly (C)';
end
for y = 1:yn
    subplot(yn,xn,1+(y-1)*xn)
    if zoom
        yyaxis left
    end
    ylabel(lbl)
end

if(tosave)
    savefig(2, ['figures/', variable, '/', scenarios{1}, '_Fig2_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N)]);
else
    fprintf(['figures/', variable, '/', scenarios{1}, '_Fig2_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '\n']);
end
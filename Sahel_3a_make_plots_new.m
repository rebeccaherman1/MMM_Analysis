historical = true;
tosave = true;
start_year = 1901;
end_year = 2013;
variable = 'pr';
realm = 'amip';%'cmip6';

switch realm
    case 'cmip6'
        scenarios = {'cmip6_h','cmip6_a','cmip6_n','cmip6_g'};%'amip'};%,; 
        scenario_names = {'ALL 6', 'AA 6', 'NAT 6', 'GHG 6'};
        pC = true;
        scenario_colors = {'b', 'm', 'r', 'g', 'c'};
        long_colors = {'blue', 'magenta', 'red', 'green', 'cyan'};
    case 'cmip5'
        scenarios = {'h', 'a', 'n', 'g'};
        scenario_names = {'ALL', 'AA', 'NAT', 'GHG'};
        pC = true;
        scenario_colors = 'bmrgc';
        long_colors = {'blue', 'magenta', 'red', 'green', 'cyan'};
    case 'amip'
        scenarios = {'cmip6_r', 'v'}; %TODO: add other amip scenarios...
        scenario_names = {'amip-hist', 'amip-piForcing'};
        pC = false;
        scenario_colors = {[0, 127, 0]/255, [1,.7,0]};
        long_colors = {'green', 'yellow'};
    otherwise
        fprintf('undefined realm!')
end
zoom = strcmp(variable, 'pr') && ~strcmp(realm, 'amip');
As = cell(length(scenarios),1);

fl = '';%'last';
N=500;

global ref_T_years; 

obs = load(['data/', variable, '/observations.mat']);

timeframe = (obs.T >= start_year & obs.T <= end_year);
ref_T_years = obs.T(timeframe); %TODO this is inconsistent
start_year = max(start_year, ref_T_years(1));
end_year = min(end_year, ref_T_years(end));
var_unstandardized = obs.var(:,timeframe,:);
var_anomaly = var_unstandardized - mean(var_unstandardized,2);
var_std = std(var_unstandardized,0,2);
var_standardized = var_anomaly./var_std;

figure(2); hold off; clf; hold on; 
if(strcmp(variable, 'ts'))
    g_GM = load(['data/', variable, '/', scenarios{contains(scenarios, 'g')}, '_GM.mat']);
    gt = single(g_GM.time);
    g_tf = ismember(gt, ref_T_years);
    no_ghg = var_anomaly - (g_GM.MMM(:,g_tf,:) - mean(g_GM.MMM(:,g_tf,:), 2));
    no_ghg(:,:,3) = var_anomaly(:,:,3);
end

N_indiv = nan(1,length(scenarios));
for j = 1:length(scenarios)
    scenario = char(scenarios(j));
    fprintf("Accessing historical scenario %s\n", scenario);
    GM = load(['data/', variable, '/', scenario, '_GM.mat']);
    A = load(['Analysis/', variable, '/', scenario, '_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.mat']);
    As{j} = A;
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
        if(strcmp(variable, 'ts'))
            r_no_ghg = permute(diag(corr(permute(mmm, [2,3,1]), permute(no_ghg, [2,3,1]))),[2,3,1]);
            e_no_ghg = mean((mmm-no_ghg).^2,2).^.5./std(no_ghg,0,2);
        end
    end
    
    %which PC do we want?
    for i=1:I
        if(I==1)
            switch realm
                case {'cmip5', 'cmip6'}
                    yn = length(scenarios)^.5; xn = yn;
                    subplot(yn, xn, j); 
                case 'amip'
                    yn = length(scenarios)+1; xn = 1;
                    subplot(yn, xn, j)
                otherwise
                    yn = 1; xn = 1;
            end
        else
            yn = length(scenarios); xn = I;
            subplot(yn,xn, I*(j-1)+i)
        end
        plot(ref_T_years, zeros(size(ref_T_years)), 'k--'); 
        hold on; set(gca,'FontSize',15); 
        if(strcmp(variable, 'ts') && (contains(scenario, 'a') || contains(scenario, 'n')) && i<3)
            p_actual = plot(ref_T_years, no_ghg(:,:,i), '-', 'color', [.5 0 .6]);%, 'LineWidth', 2);
            r_ttl = r_no_ghg(:,:,i); e_ttl = e_no_ghg(:,:,i);
            fprintf('no ghg!')
        else
            p_actual = plot(ref_T_years, var_anomaly(:,:,i), 'k-');%, 'LineWidth', 2);
            r_ttl = r(:,:,i); e_ttl = rmsd(:,:,i);
        end
        if(I==1)
            title([scenario_names{j}, ', N=', num2str(N_indiv(j)), ', r=', num2str(r_ttl, '%5.2f'), ', rmse=', num2str(rmsd, '%5.2f')], 'color', scenario_colors{j})
            ylabel('Precipitation Anomaly (mm/day)')
        else
            if(i==1)
                ylabel([scenario_names{j}, ', N=', num2str(N_indiv(j))], 'color', scenario_colors{j})
            end
            ttl = ['\color{', long_colors{j}, '} r=', num2str(r_ttl, '%5.2f'), ', rmse=', num2str(e_ttl, '%5.2f')];
            if(j==1)
                ttl = [ttl, ' \color{black}', A.indices{i}];
            end
            title(ttl)
        end 
        if(zoom)
            yyaxis right
            set(gca, 'ycolor', scenario_colors{j})
            ylim([-.375,.375])
        else
            %set(gca, 'ycolor', scenario_colors{j})
            ylim([-1.5,1.5])
        end
        %pC
        if(pC)
            p_stds = fill([ref_T_years';flipud(ref_T_years')],[pC_down(:,:,i)';flipud(pC_up(:,:,i)')],'y','FaceAlpha', .3 ,'linestyle','none');
            p_stds.HandleVisibility = 'off';
        end
        if(historical)
            hold on;
            p_stds = fill([ref_T_years';flipud(ref_T_years')],[down(:,:,i)';flipud(up(:,:,i)')],scenario_colors{j},'FaceAlpha', .3 ,'linestyle','none');
            p_stds.HandleVisibility = 'off';
            p_mmm = plot(ref_T_years,mmm(:,:,i),'-', 'Color', scenario_colors{j});%, 'LineWidth', 2);%, 'DisplayName', );
            xlim([start_year, end_year]);
        end
        %still need this for pr... and it has to be at the end, it would
        %seem.
    end
end

if(strcmp(realm, 'amip'))
    R.diffs = As{1}.historical_bootstrapped.b_means - As{2}.historical_bootstrapped.b_means;
    %R.MMM = mean(R.diffs, 1);
    R.MMM = As{1}.MMM.MMM - As{2}.MMM.MMM; 
    R.MMM = R.MMM - mean(R.MMM, 2);
    [R.low, R.high] = confidence_interval(R.diffs);
    subplot(yn, xn, yn); hold on; set(gca,'FontSize',15); 
    plot(ref_T_years, zeros(size(ref_T_years)), 'k--'); 
    r_stds = fill([ref_T_years';flipud(ref_T_years')],[R.low(:,:,i)';flipud(R.high(:,:,i)')],[191, 0, 191]/255, 'FaceAlpha', .3 ,'linestyle','none');
    r_stds.HandleVisibility = 'off';
    p_mmm = plot(ref_T_years,R.MMM(:,:,i),'-', 'Color', [126, 47, 142]/255);%, 'LineWidth', 2);%, 'DisplayName', );
    xlim([start_year, end_year]); ylim([-1,1])
    n = load('data/pr/cmip6_n_GM.mat');
    mmm = n.MMM - mean(n.MMM);
    plot(ref_T_years, mmm(ismember(single(n.time(1,:)), ref_T_years)), 'r-')
    title('Implied Fast Component', 'Color', [126, 47, 142]/255)

end

if(tosave)
    savefig(2, ['figures/', variable, '/', realm, '_Fig2_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N)]);
    saveas(2, ['figures/', variable, '/', realm, '_Fig2_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N)], 'png');
else
    fprintf(['figures/', variable, '/', realm, '_Fig2_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '\n']);
end
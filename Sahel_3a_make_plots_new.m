clear
historical = true;
tosave = true;
start_year = 1901;
anomaly_years = 1901:1950;
variable = 'ts';
realm = 'cmip5';
%TODO add NARI for amip figures.

%TODO: (automate adding cmip5.)
%maybe, add ALL MMM to IFC.

switch realm
    case 'cmip6'
        scenarios = {'cmip6_h','cmip6_a','cmip6_n','cmip6_g'};%'amip'};%,; 
        scenario_names = {'ALL 6', 'AA 6', 'NAT 6', 'GHG 6'};
        pC = true;
        scenario_colors = {'b', 'm', [0.60,0.20,0.00], [0.00,0.80,0.00], 'c'};
        long_colors = {'blue', 'magenta', 'red', 'green', 'cyan'};
        end_year = 2014;
    case 'cmip5'
        scenarios = {'h', 'a', 'n', 'g'};
        scenario_names = {'ALL', 'AA', 'NAT', 'GHG'};
        pC = true;
        scenario_colors = {[0.00,0.45,0.74], [0.75,0.00,0.75], 'r', [0.47,0.67,0.19], 'c'};
        long_colors = {'blue', 'magenta', 'red', 'green', 'cyan'};
        end_year = 2003;
    case 'amip'
        scenarios = {'amip-piF', 'amip-hist', 'cmip6_fast'}; %TODO: add other amip scenarios...
        scenario_names = {'amip-piForcing','amip-hist', 'Implied Fast Component'};
        pC = false;
        scenario_colors = {[1,.7,0],[0, 127, 0]/255, [126, 47, 142]/255};
        end_year = 2014;
    otherwise
        fprintf('undefined realm!')
end
zoom = strcmp(variable, 'pr') && ~strcmp(realm, 'amip'); %I put a special case for 'cmip6_fast' scenario at the end
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
var_anomaly = var_unstandardized - mean(var_unstandardized(:,ismember(ref_T_years, anomaly_years),:),2);
var_std = std(var_unstandardized,0,2);
var_standardized = var_anomaly./var_std;

if(~strcmp(realm, 'amip'))
    gA = load(['Analysis/', variable, '/', scenarios{contains(scenarios, 'g')}, '_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.mat']);
    gt = start_year:end_year;
    g_tf = ismember(gt, ref_T_years);
    a_g_tf = ismember(gt, anomaly_years);
    GHG_sub = (gA.MMM.MMM(:,g_tf,:) - mean(gA.MMM.MMM(:,a_g_tf,:), 2));
    no_ghg = var_anomaly - GHG_sub;
    if(strcmp(variable, 'ts'))
        no_ghg(:,:,3) = var_anomaly(:,:,3);
    end
end

%% Make Fig2

figure(2); hold off; clf; hold on; 

N_indiv = nan(1,length(scenarios));
for j = 1:length(scenarios)
    scenario = char(scenarios(j));
    fprintf("Accessing historical scenario %s\n", scenario);
    A = load(['Analysis/', variable, '/', scenario, '_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.mat']);
    if(strcmp(variable, 'ts'))
        A.indices = {'NA', 'GT', 'NARI'};
    end
    As{j} = A;
    if(isfield(A, 'indiv'))
        N_indiv(j) = length(A.indiv.models);
    else
        N_indiv(j) = 0;
    end
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
        mmm_v = A.MMM.MMM; mmm_anomaly = mmm_v-mean(mmm_v,2); r = A.MMM.r; rmsd = A.MMM.e;
        anomaly_diff = mean(mmm_anomaly(:,ismember(ref_T_years, anomaly_years),:),2);
        mmm = mmm_anomaly-anomaly_diff;
        %the historical_bootstrapped are actual anomalies, hence the need
        %for "anomaly_diff" as defined.
        up   = A.historical_bootstrapped.high - anomaly_diff; 
        down = A.historical_bootstrapped.low - anomaly_diff;  
        if(strcmp(variable, 'ts'))
            r_no_ghg = permute(diag(corr(permute(mmm_anomaly, [2,3,1]), permute(no_ghg, [2,3,1]))),[2,3,1]);
            e_no_ghg = mean((mmm_anomaly-no_ghg).^2,2).^.5./std(no_ghg,0,2);
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
                    yn = length(scenarios); xn = 1;
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
        
        %formatting stuffs
        if(I==1)
            title([scenario_names{j}, ', N=', num2str(N_indiv(j)), ', r=', num2str(r_ttl, '%5.2f'), ', sRMSE=', num2str(rmsd, '%5.2f')], 'color', scenario_colors{j})
            ylabel('Precipitation Anomaly (mm/day)')
        else
            if(i==1)
                ylabel([scenario_names{j}, ', N=', num2str(N_indiv(j))], 'color', scenario_colors{j})
            end
            ttl = ['r=', num2str(r_ttl, '%5.2f'), ', sRMSE=', num2str(e_ttl, '%5.2f')];%\color{', long_colors{j}, '} 
            if(j==1)
                ttl = [ttl, ' \color{black}', A.indices{i}];
            end
            title(ttl, 'Color', scenario_colors{j})
        end 
        %left ordinates
        if(zoom)
            yl = 1.25;
            ylim([-yl, yl])
        elseif(~strcmp(variable,'ts'))
            %set(gca, 'ycolor', scenario_colors{j})
            ylim([-1.5,1.5])
        end
        %right ordinates
        if(zoom || strcmp(scenario, 'cmip6_fast'))% || (I~=1 && strcmp(A.indices(i), 'NARI'))
            yyaxis right
            set(gca, 'ycolor', scenario_colors{j})
        end
        if(strcmp(scenario, 'cmip6_fast'))
            ylim([-1,1]);
        elseif(zoom || I~=1 && strcmp(A.indices(i), 'NARI'))
            ylim([-.5, .5])
        elseif(strcmp(variable, 'ts'))
            %ylim([-.85,1.15]);
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
        %{
        if(strcmp(variable, 'ts'))
            if i==3
                yyaxis left
            end
            yl = ylim;
            ylim([floor(yl(1)), ceil(yl(2))])
        end
        %}
    end
end

if(strcmp(realm, 'amip'))
    NAT = load(['Analysis/pr/cmip6_n_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.mat']);
    mmm = NAT.MMM.MMM;
    mmm = mmm - mean(mmm(:,ismember(ref_T_years, anomaly_years),:),2);
    plot(ref_T_years, mmm, 'r-')
    corr(mmm', A.MMM.MMM')
    
    ALL = load(['Analysis/pr/cmip6_h_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.mat']);
    mmm = ALL.MMM.MMM;
    mmm = mmm - mean(mmm(:,ismember(ref_T_years, anomaly_years),:),2);
    plot(ref_T_years, mmm, 'b-')
    corr(mmm', A.MMM.MMM')
    
    sst_obs = load('data/ts/observations.mat'); NARI = sst_obs.var(:,:,3);
    NARI_color = [0.07,0.62,1.00];
    subplot(yn, xn, 1); yyaxis right; set(gca, 'ycolor', NARI_color); ylabel('NARI (C)');
    plot(sst_obs.T, NARI, '-', 'Color', NARI_color, 'DisplayName', 'NARI')
end
%}
if(tosave)
    savefig(2, ['figures/', variable, '/', realm, '_Fig2_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N)]);
    saveas(2, ['figures/', variable, '/', realm, '_Fig2_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N)], 'png');
else
    fprintf(['figures/', variable, '/', realm, '_Fig2_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '\n']);
end


%% Make Booth Figure
% This figure is not saved.
if(strcmp(variable, 'ts'))
    GHG_sub = GHG_sub(:,:,1);
    lin_trend = GHG_sub - detrend(GHG_sub);
    %not as dramatic without the extra models in cmip5

    figure(6); clf;
    plot(ref_T_years, no_ghg(:,:,1), 'k-', 'DisplayName', 'Observations - GHG MMM'); hold on;
    plot(ref_T_years, var_anomaly(:,:,1) - lin_trend, 'k-.', 'DisplayName', 'Observations - Linear Trend');

    A = load(['analysis/', variable, '/', scenarios{1}, '_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.mat']);
    mmm_v = A.MMM.MMM(:,:,1); mmm_anomaly = mmm_v-mean(mmm_v,2); 
    anomaly_diff = mean(mmm_anomaly(:,ismember(ref_T_years, anomaly_years),:));
    mmm = mmm_anomaly-anomaly_diff;
    up   = A.historical_bootstrapped.high(:,:,1) - anomaly_diff; 
    down = A.historical_bootstrapped.low(:,:,1) - anomaly_diff;  

    mmm_noghg = mmm - GHG_sub;
    offset = mean(mmm_noghg(:,a_g_tf,:),2);
    mmm_noghg = mmm_noghg - offset;
    up_noghg = up - GHG_sub - offset;
    down_noghg = down - GHG_sub - offset;

    p_stds = fill([ref_T_years';flipud(ref_T_years')],[down_noghg(:,:,1)';flipud(up_noghg(:,:,1)')],'b','FaceAlpha', .3 ,'linestyle','none');
    hold on; p_stds.HandleVisibility = 'off';
    plot(ref_T_years,mmm_noghg(:,:,1),'-', 'Color', 'b', 'DisplayName', 'ALL - GHG MMM');

    mmm = mmm - lin_trend;
    offset = mean(mmm(:,a_g_tf,:),2);
    mmm = mmm - offset;
    up = up - lin_trend - offset;
    down = down - lin_trend - offset;

    p_stds = fill([ref_T_years';flipud(ref_T_years')],[down(:,:,1)';flipud(up(:,:,1)')],[255,153,0]/255,'FaceAlpha', .3 ,'linestyle','none');
    p_stds.HandleVisibility = 'off';
    plot(ref_T_years,mmm(:,:,1),'-', 'Color', [255,153,0]/255, 'DisplayName', 'ALL - Linear Trend');
    
    %AA
    A = load(['analysis/', variable, '/', scenarios{2}, '_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.mat']);
    mmm_v = A.MMM.MMM(:,:,1); mmm_anomaly = mmm_v-mean(mmm_v,2); 
    anomaly_diff = mean(mmm_anomaly(:,ismember(ref_T_years, anomaly_years),:));
    mmm = mmm_anomaly-anomaly_diff;
    up   = A.historical_bootstrapped.high(:,:,1) - anomaly_diff; 
    down = A.historical_bootstrapped.low(:,:,1) - anomaly_diff;  
    
    p_stds = fill([ref_T_years';flipud(ref_T_years')],[down(:,:,1)';flipud(up(:,:,1)')],'m','FaceAlpha', .3 ,'linestyle','none');
    p_stds.HandleVisibility = 'off';
    plot(ref_T_years,mmm(:,:,1),'-', 'Color', 'm', 'DisplayName', 'AA');

    yl = ylim; YL = mean(abs(yl));
    ylim([-YL, YL]);
    xlim([start_year, end_year])

    legend('show')
end
clear
historical = true;
tosave = false;
start_year = 1901;
anomaly_years = 1901:1920;
variable = 'ts';
realm = 'cmip6';
start_month = 7;
end_month = 9;
%TODO add NARI for amip figures.

%TODO: (automate adding cmip5.)
%maybe, add ALL MMM to IFC.

%TODO make new case for CMIP6 (and CMIP5?) fast.
fltr = 20;
gry = [0.65,0.65,0.65]; ls = '-'; lw = 2;
v_xticks = [1912,1920,1940,1963,1982,1991,2000];
v_xticklabels = {'- Katmai (1912)','1920','1940','+ Agung (1963)',...
    '- El Chichon (1982)','+ Pinatubo (1991)','2000'};
v_xta = 40;

switch realm
    case 'cmip6'
        %TODO this isn't generally for fast... make it its own special
        %case.
        %scenarios = {'cmip6_hfast','cmip6_afast','cmip6_nfast','cmip6_gfast'};%'amip'};%,; 
        scenarios = {'cmip6_h','cmip6_a','cmip6_n','cmip6_g'};%'amip'};%,; 
        scenario_names = {'a. ALL', 'b. AA', 'c. NAT', 'd. GHG'};
        pC = true;
        scenario_colors = {'b', 'm', [0.60,0.20,0.00], [0.00,0.80,0.00], 'c'};
        long_colors = {'blue', 'magenta', 'red', 'green', 'cyan'};
        end_year = 2014;
    case 'cmip5'
        %scenarios = {'hfast', 'afast', 'nfast', 'gfast'};
        scenarios = {'h', 'a', 'n', 'g'};
        scenario_names = {'ALL', 'AA', 'NAT', 'GHG'};
        pC = true;
        scenario_colors = {[0.00,0.45,0.74], [0.75,0.00,0.75], 'r', [0.47,0.67,0.19], 'c'};
        long_colors = {'blue', 'magenta', 'red', 'green', 'cyan'};
        end_year = 2003;
        ls = ':'; lw = 2;
    case 'amip'
        scenarios = {'amip-piF', 'amip-hist', 'cmip6_fast'}; %TODO: add other amip scenarios...
        scenario_names = {'amip-piForcing','amip-hist', 'Implied Fast Component'};
        pC = false;
        scenario_colors = {[1.00,0.54,0.16],[0, 127, 0]/255, [126, 47, 142]/255};
        end_year = 2014;
        anomaly_years = 1901:1920;
    otherwise
        fprintf('undefined realm!')
end
zoom = 0;%strcmp(variable, 'pr') && ~strcmp(realm, 'amip') && ~any(contains(scenarios, 'fast')); %I put a special case for 'cmip6_fast' scenario at the end
As = cell(length(scenarios),1);

fl = '';%'last';
N=500;

global ref_T_years; 

%{
F = load('Analysis/pr/cmip6_fast_1901-2014_N500.mat');
F.T = 1901:2014;
ref_T_years = start_year:end_year;
var_unstandardized = F.MMM.MMM(ismember(F.T,ref_T_years));
var_smthd = smooth(var_unstandardized,fltr)';
%}
O = load(sprintf('data/%s/7-9/observations.mat', variable));
ref_T_years = O.T(ismember(O.T, start_year:end_year));

%SMOOTH IT UP HERE
var_smthd = smoothdata(O.var, 2, 'movmean', fltr);
var_smthd = var_smthd(:, ismember(O.T, start_year:end_year),:);
var_unstandardized = O.var(:,ismember(O.T, start_year:end_year),:);
%}
%obs = load(make_data_filename(variable, start_month, end_month, 'observations'));
%timeframe = (obs.T >= start_year & obs.T <= end_year);
%ref_T_years = obs.T(timeframe); %TODO this is inconsistent
start_year = max(start_year, ref_T_years(1));
end_year = min(end_year, ref_T_years(end));
%var_unstandardized = obs.var(:,timeframe,:);
var_anomaly = var_unstandardized - mean(var_unstandardized(:,ismember(ref_T_years, anomaly_years),:),2);
smthd_anomaly = var_smthd - mean(var_smthd(:, ismember(ref_T_years, anomaly_years),:),2);
var_std = std(var_unstandardized,0,2);
var_standardized = var_anomaly./var_std;

if(~strcmp(realm, 'amip') && ~contains(scenarios{1}, 'fast'))
    gA = load(make_analysis_filename(variable,scenarios{contains(scenarios, 'g')}, start_year, end_year, N));
    gt = start_year:end_year;
    g_tf = ismember(gt, ref_T_years);
    a_g_tf = ismember(gt, anomaly_years);
    %always present smoothed versions
    GHG_sub = (gA.MMM.MMM(:,g_tf,:) - mean(gA.MMM.MMM(:,a_g_tf,:), 2));
    no_ghg = var_anomaly - GHG_sub;
    if(strcmp(variable, 'ts'))
        no_ghg(:,:,3) = var_anomaly(:,:,3);
    end
end

%var_anomaly and no_ghg are NOT SMOOTHED. must do later.

%% Make Fig2

figure(2); hold off; clf; hold on; 

N_indiv = nan(1,length(scenarios));
for j = 1:length(scenarios)
    scenario = char(scenarios(j));
    smth = strcmp(scenario, 'a') || strcmp(scenario, 'afast') || contains(scenario, {'_a', 'g', 'amip', '_f'});%, 'fast'});
    fprintf("Accessing historical scenario %s\n", scenario);
    A = load(make_analysis_filename(variable, scenario, start_year, end_year, N));
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

    DU = {'low', 'high'};
    MMM = 'MMM';
    if(smth)
        DU = arrayfun(@(X) append(X, '_s'), DU);
        MMM = append(MMM, '_s');
    end
    if(pC)
        if(strcmp(fl, 'last'))
            pC_down = A.piC_last.low;    
            pC_up = A.piC_last.high;     
            pC_scale = mean(std(A.piC_last.r_means, 0, 2));         
        else
            pC_down = nanmean(A.piC_resampled_bootstrapped.(DU{1}));
            pC_up   = nanmean(A.piC_resampled_bootstrapped.(DU{2}));
            pC_scale = nanmean(std(A.piC_resampled_bootstrapped.r_means, 0, 2));
        end
    end
    if(historical)
        mmm_v = A.(MMM).MMM; mmm_anomaly = mmm_v-mean(mmm_v,2); 
        %if(not(strcmp(realm, 'amip')))
            r_nm = ', r_L_F=';
            sRMSE_nm = ', sRMSE_L_F=';
            r = A.MMM_s.r; rmsd = A.MMM_s.e;
        %else
        %    r_nm = ', r=';
        %    sRMSE_nm = ', sRMSE=';
        %    r = A.MMM.r; rmsd = A.MMM.e;
        %end
        if(contains(scenario, 'n'))
            anomaly_diff = mean(mmm_anomaly(:,20:60,:),2);
        else
            anomaly_diff = mean(mmm_anomaly(:,ismember(ref_T_years, anomaly_years),:),2);
        end
        mmm = mmm_anomaly-anomaly_diff;
        %the historical_bootstrapped are actual anomalies, hence the need
        %for "anomaly_diff" as defined.
        up   = A.historical_bootstrapped.(DU{2}) - anomaly_diff; 
        down = A.historical_bootstrapped.(DU{1}) - anomaly_diff; 
        if(smth)
            mmm = smoothdata(mmm, 2, 'movmean', 5);
            up = smoothdata(up, 2, 'movmean', 5);
            down = smoothdata(down, 2, 'movmean', 5);
        end
        if(strcmp(variable, 'ts'))
            if(smth)
                no_ghg_use = smoothdata(no_ghg, 2, 'movmean', fltr);
            else
                no_ghg_use = no_ghg;
            end
            r_no_ghg = permute(diag(corr(permute(mmm_anomaly, [2,3,1]), permute(no_ghg_use, [2,3,1]))),[2,3,1]);
            e_no_ghg = mean((mmm_anomaly-no_ghg_use).^2,2).^.5./std(no_ghg_use,0,2);
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
        smthd = smooth(smthd_anomaly(:,:,i), 5);
        tc = 'k'; r_ttl = r(:,:,i); e_ttl = rmsd(:,:,i);
        if(strcmp(variable, 'ts') && contains(scenario, 'a') && i<3)
            tp = no_ghg_use(:,:,i)';
            tc = [1.00,0.55,0.31];
            r_ttl = r_no_ghg(:,:,i); e_ttl = e_no_ghg(:,:,i);
            fprintf('no ghg!')
        elseif(smth)
            tp = smthd_anomaly(:,:,i)';
            %just HF var
        elseif(contains(scenario, 'n'))
            tp = var_anomaly(:,:,i)-smthd_anomaly(:,:,i);
            if(strcmp(variable, 'ts'))
                r_ttl = r_no_ghg(:,:,i); e_ttl = e_no_ghg(:,:,i);
            end
            tc = gry;
        else
            tp = var_anomaly(:,:,i);
        end
        if(smth)
            tp = smoothdata(tp, 1, 'movmean', 5);
        end
        if(strcmp(scenario, 'h') || contains(scenario, '_h') || strcmp(scenario, 'hfast'))
            plot(ref_T_years, smthd, '-', 'color', 'k');%*.8/.65
            tc = gry;
        end
        p_actual_s = plot(ref_T_years, tp, '-', 'color', tc);%, 'LineWidth', lw);
        %formatting stuffs
        if(I==1)
            title([scenario_names{j}, r_nm, num2str(r_ttl, '%5.2f'), sRMSE_nm, num2str(rmsd, '%5.2f')], 'color', scenario_colors{j})%', N=', num2str(N_indiv(j)),
            ylabel('Precipitation Anomaly (mm/day)')
        else
            if(i==1)
                ylabel(scenario_names{j}, 'color', scenario_colors{j})%', N=', num2str(N_indiv(j))
            end
            ttl = ['r=', num2str(r_ttl, '%5.2f'), sRMSE_nm, num2str(e_ttl, '%5.2f')];%\color{', long_colors{j}, '} 
            if(j==1)
                ttl = {['{\color{black}', A.indices{i},'}'], ttl};
            end
            title(ttl, 'Color', scenario_colors{j})
        end 
        N_indiv(j)
        if(contains(scenario, {'n'}) && ~strcmp(variable, 'ts'))%, 'fast'}))
            xticks(v_xticks);
            xticklabels(v_xticklabels);
            xtickangle(v_xta);
        end
        %left ordinates
        if(zoom)
            yl = 1.25;
            ylim([-yl, yl])
            yyaxis right
            set(gca, 'ycolor', scenario_colors{j})
        elseif(~strcmp(variable,'ts') && ~any(contains(scenarios, 'fast')))
            %set(gca, 'ycolor', scenario_colors{j})
            %ylim([-1.5,1.5])
            ylim([-.7,.7])
        end

        if(strcmp(realm, 'amip'))
            if(smth)
                ylim([-0.5000,0.5849])
            else
                ylim([-1.5,1.5]);
            end
        end
        
        if(zoom || I~=1 && strcmp(A.indices(i), 'NARI'))
            ylim([-.5, .5])
        elseif(strcmp(variable, 'ts'))
            %ylim([-.85,1.15]);
        end

        %pC
        if(pC)
            %p_stds = fill([ref_T_years';flipud(ref_T_years')],[pC_down(:,:,i)';flipud(pC_up(:,:,i)')],'y','FaceAlpha', .3 ,'linestyle','none');
            %p_stds.HandleVisibility = 'off';
            plot([ref_T_years(1), 2014], [1,1]*pC_up(:,:,i), 'k', 'LineStyle', ls, 'HandleVisibility', 'off')
            plot([ref_T_years(1), 2014], [1,1]*pC_down(:,:,i), 'k', 'LineStyle', ls, 'HandleVisibility', 'off')
        end
        if(historical)
            %I've already smoothed it now!
            %{
            if(smth)
                mmm(:,:,i) = smooth(mmm(:,:,i), fltr);
                down(:,:,i) = smooth(down(:,:,i), fltr);
                up(:,:,i) = smooth(up(:,:,i), fltr); 
            end
            %}
            hold on;
            p_stds = fill([ref_T_years';flipud(ref_T_years')],[down(:,:,i)';flipud(up(:,:,i)')],scenario_colors{j},'FaceAlpha', .3 ,'linestyle','none');
            p_stds.HandleVisibility = 'off';
            p_mmm = plot(ref_T_years,mmm(:,:,i),'LineStyle', ls, 'Color', scenario_colors{j}, 'LineWidth', lw);%, 'DisplayName', );
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
    %{
    NAT = load(make_analysis_filename('pr','cmip6_n', start_year, end_year, N));
    mmm = NAT.MMM.MMM;
    mmm = mmm - mean(mmm(:,ismember(ref_T_years, anomaly_years),:),2);
    plot(ref_T_years, mmm, 'r-')
    corr(mmm', A.MMM.MMM')
    %}
    %{
    ALL = load(make_analysis_filename('pr', 'cmip6_h', start_year, end_year, N));
    if(smth)
        mmm = ALL.MMM_s.MMM;
    else
        mmm = ALL.MMM.MMM;
    end
    mmm = mmm - mean(mmm(:,ismember(ref_T_years, anomaly_years),:),2);
    plot(ref_T_years, mmm, 'b-')
    corr(mmm', A.MMM.MMM')
    %}
    sst_obs = load(make_data_filename('ts', start_month, end_month, 'observations')); 
    if(0)
        NARI = smooth(smooth(sst_obs.var(:,:,3), fltr),5);
    else
        NARI = sst_obs.var(:,:,3);
    end
    NARI = NARI(ismember(sst_obs.T, ref_T_years));
    NARI = NARI - mean(NARI(ismember(ref_T_years, anomaly_years)));
    NARI_color = [0.30,0.75,0.93];%0.07,0.62,1.00];
    subplot(yn, xn, 1); yl = ylim;
    yyaxis right; set(gca, 'ycolor', NARI_color); ylabel('NARI Anomaly (C)');
    ylim(yl./.87) %mean(NARI(ismember(ref_T_years, anomaly_years)))
    plot(ref_T_years, NARI, '-', 'Color', NARI_color, 'DisplayName', 'NARI')
end
%}
if(tosave)
    savefig(2, ['figures/', variable, '/LF/', realm, '_Fig2_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N)]);
    saveas(2, ['figures/', variable, '/LF/', realm, '_Fig2_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N)], 'png');
else
    fprintf(['figures/', variable, '/LF/', realm, '_Fig2_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '\n']);
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

    A = load(make_analysis_filename(variable, scenarios{1}, start_year, end_year, N));
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
    A = load(make_analysis_filename(variable, scenarios{2}, start_year, end_year, N));
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

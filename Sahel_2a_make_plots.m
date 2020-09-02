save = true;
variable = 'pr';
start_year=1901;
end_year=2014;
%should work with TS if I pick a basin now, but no need to make this figure
%for that.

scenarios = {'v'};%'cmip6_h'};%'r'};%'a6'};%'amip'};%'h'};%,"historicalAerosol","historicalNat","historicalGHG"  "historicalMisc", 

global ref_T_years

obs = load(['data/', variable, '/observations.mat']);
%this is still short
%cru = load(['data/', variable, '/CRU_data.mat']);%ncread('data/Jul-Sep/CRU_data.nc', 'aprod'); %mm/month *month/day
%cru=cru.prcp;
timeframe_obs = (obs.T >= start_year & obs.T <= end_year);
ref_T_years = obs.T(timeframe_obs);
prcp = obs.var(:,ismember(obs.T, ref_T_years),:); %cru=cru(timeframe_obs);
prcp_anomaly = prcp - mean(prcp); %cru_anomaly = cru - mean(cru);
prcp_standardized = prcp_anomaly./std(prcp,0,2); %cru_standardized = cru_anomaly/std(cru);
    
for i = 1:length(scenarios)
    scenario = scenarios{i};
    fprintf("Accessing historical scenario %s\n", scenario);
    G = load(['data/', variable, '/',scenario,'_GM.mat']);
    timeframe_m = ismember(single(G.time), ref_T_years);
    
    MMM = G.MMM(:,timeframe_m,:);                
    MMM_anomaly = MMM - mean(MMM,2); 
    MMM_standardized = MMM_anomaly./std(MMM_anomaly,0,2);

    GM = G.GMs(:,timeframe_m,:);
    GM_anomalies = GM - mean(GM,2);   
    GM_standardized = GM_anomalies./std(GM_anomalies, 0, 2);
    
    I = load(['data/', variable, '/',scenario,'_all.mat']);
    runs = I.runs(:,timeframe_m); 
    runs_anomalies = runs - mean(runs,2);   
    runs_standardized = runs_anomalies./std(runs_anomalies, 0, 2);

    figure(1); hold off; clf;

    subplot(2,2,1); hold on; %yikes 
    set(gca,'FontSize',16); %set(gca,'LineWidth',2);
    title('a. Standardized');
    ylabel("Standardized Precipitation Anomaly");
    plot(ref_T_years,runs_standardized, '-', 'Color', [.8, .8, 1], 'LineWidth', .1, 'HandleVisibility', 'off');
    plot(ref_T_years,GM_standardized, 'c-', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(ref_T_years,MMM_standardized,'b-', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(ref_T_years,prcp_standardized, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    %plot(ref_T_years,cru_standardized, '--', 'Color', .8*[1,1,1], 'HandleVisibility', 'off');
    ylim([-4,4]); xlim([start_year, end_year]);

    subplot(2,2,2); hold on; %yikes again
    set(gca,'FontSize',16); %set(gca,'LineWidth',2);
    title('b. Anomalies');
    ylabel("Precipitation Anomaly (mm/day)"); 
    p_runs_s = plot(ref_T_years,runs_anomalies, '-', 'Color', [.8, .8, 1], 'LineWidth', .1);
    p_gm_s = plot(ref_T_years,GM_anomalies, 'c-', 'LineWidth', 1);
    p_mmm_s = plot(ref_T_years,MMM_anomaly,'b-', 'LineWidth', 2);
    p_actual_s = plot(ref_T_years,prcp_anomaly, 'k-', 'LineWidth', 2);
    %p_actual_s_cru = plot(ref_T_years,cru_anomaly, '--', 'Color', .8*[1,1,1]);
    xlim([start_year, end_year]); ylim([-2,2]);

    legend([p_runs_s(1), p_gm_s(1), p_mmm_s, p_actual_s], 'Runs', 'IMs', 'MMM', 'GPCC', 'Location', 'south');%, p_actual_s_cru]'northwest'); 'CRU', %

    if(save)
        savefig(['figures/', variable, '/', scenario, '_Fig1p_', num2str(ref_T_years(1)), '-', num2str(ref_T_years(end))]);
    end
end

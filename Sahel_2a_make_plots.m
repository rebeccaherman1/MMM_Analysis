save = true;
short=false;%true;

scenarios = {'v'};%'r'};%'a6'};%'amip'};%'h'};%,"historicalAerosol","historicalNat","historicalGHG"  "historicalMisc", 
if(short || strcmp(scenarios, 'amip'))
    start_year = 1950;% 
else
    start_year = 1901;%
end
end_year=2003;

global ref_T_years ref_T_amip;

obs = load('model_output/historical_precipitation.mat');
cru = load('model_output/CRU_data.mat');%ncread('model_output/Jul-Sep/CRU_data.nc', 'aprod'); %mm/month *month/day
cru=cru.prcp;
%cru = cru/mean([31,31,30]); now we did this before we saved it
timeframe_obs = (obs.T >= start_year & obs.T <= end_year);
ref_T_amip = start_year:end_year;
ref_T_years = obs.T(timeframe_obs);
prcp = obs.prcp(timeframe_obs); cru=cru(timeframe_obs);
prcp_anomaly = prcp - mean(prcp); cru_anomaly = cru - mean(cru);
prcp_standardized = prcp_anomaly/std(prcp); cru_standardized = cru_anomaly/std(cru);
    
for i = 1:length(scenarios)
    scenario = char(scenarios(i));
    fprintf("Accessing historical scenario %s\n", scenario);

    if(short || strcmp(scenario, 'amip'))
        timeframe_m = ref_T_amip>=start_year & ref_T_amip<=end_year;
    else
        timeframe_m = timeframe_obs;
    end
    
    G = load(['model_output/',scenario,'_GM.mat']);
    MMM = G.MMM(timeframe_m); %MMM = MMM(1:99);                
    MMM_anomaly = MMM - mean(MMM); 
    MMM_standardized = MMM_anomaly/std(MMM_anomaly);

    GM = G.GMs(:,timeframe_m);%(:,1:99); 
    GM_anomalies = GM - mean(GM,2);   
    GM_standardized = GM_anomalies./std(GM_anomalies, 1, 2);

    I = load(['model_output/',scenario,'_all.mat']);
    runs = I.runs(:,timeframe_m);%(:,1:99); 
    runs_anomalies = runs - mean(runs,2);   
    runs_standardized = runs_anomalies./std(runs_anomalies, 1, 2);

    figure(1); hold off; clf;

    subplot(2,2,1); hold on; %yikes 
    set(gca,'FontSize',16); %set(gca,'LineWidth',2);
    title('a. Standardized');
    ylabel("Standardized Precipitation Anomaly");
    plot(ref_T_years,runs_standardized, '-', 'Color', [.8, .8, 1], 'LineWidth', .1, 'HandleVisibility', 'off');
    plot(ref_T_years,GM_standardized, 'c-', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(ref_T_years,MMM_standardized,'b-', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(ref_T_years,prcp_standardized, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(ref_T_years,cru_standardized, '--', 'Color', .8*[1,1,1], 'HandleVisibility', 'off');
    ylim([-4,4]); xlim([start_year, end_year]);

    subplot(2,2,2); hold on; %yikes again
    set(gca,'FontSize',16); %set(gca,'LineWidth',2);
    title('b. Anomalies');
    ylabel("Precipitation Anomaly (mm/day)"); 
    p_runs_s = plot(ref_T_years,runs_anomalies, '-', 'Color', [.8, .8, 1], 'LineWidth', .1);
    p_gm_s = plot(ref_T_years,GM_anomalies, 'c-', 'LineWidth', 1);
    p_mmm_s = plot(ref_T_years,MMM_anomaly,'b-', 'LineWidth', 2);
    p_actual_s = plot(ref_T_years,prcp_anomaly, 'k-', 'LineWidth', 2);
    p_actual_s_cru = plot(ref_T_years,cru_anomaly, '--', 'Color', .8*[1,1,1]);
    xlim([start_year, end_year]); ylim([-2,2]);

    legend([p_runs_s(1), p_gm_s(1), p_mmm_s, p_actual_s, p_actual_s_cru], 'Runs', 'IMs', 'MMM', 'GPCC', 'CRU', 'Location', 'south');%'northwest'); %

    if(save)
        savefig([scenario,'_Fig1']);
    end
end

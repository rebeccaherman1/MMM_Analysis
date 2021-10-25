save = false;
start_year=1901;
anomaly_years = 1901:1950;
%should work with TS if I pick a basin now, but no need to make this figure
%for that.
start_month = 7;%7
end_month = 9;%9

%TODO put colors other than blue for the different scenarios!
scenarios = {'cmip6_hfast', 'cmip6_afast', 'cmip6_nfast', 'cmip6_gfast'};%,'h','amip-hist','amip-piF'};

global ref_T_years
    
for i = 1:length(scenarios)
    scenario = scenarios{i};
    switch scenario
        case {'cmip6_h','h',}
            variables = {'pr'};%, 'ts'}; TODO not currently implemented for ts
            %for ts, we don't need the standardized anomalies. we can just
            %do anomalies, correlation, and rmse. would have to make
            %changes in the following file as well.
            clr = 'b'; mdgnd = 'c'; bckgnd = [.8, .8, 1];
        case 'amip-hist'
            variables = {'pr'};
            clr = [0, 127, 0]/255;
            mdgnd = max(min(clr*2, [1,1,1]), [.4,.4,.6]);
            bckgnd = max(clr/max(clr), [.9,.8,.8]);
            end_year = 2014;
        case 'amip-piF'
            variables = {'pr'};
            clr = [1,.7,0];
            mdgnd = max(min(clr*1.5, [1,1,1]), [.4,.6,.4]);
            bckgnd = max(clr/max(clr), [.9,.8,.8]);
    end
    if(strcmp(scenario, 'h'))
        end_year = 2003;
    else
        end_year = 2014;
    end
    fprintf("Accessing historical scenario %s\n", scenario);
    
    for v = 1:length(variables)
        variable = variables{v};
	obs = load(make_data_filename(variable, start_month, end_month, 'observations'));
        %this is still short
        %cru = load(['data/', variable, '/CRU_data.mat']);%ncread('data/Jul-Sep/CRU_data.nc', 'aprod'); %mm/month *month/day
        %cru=cru.prcp;
        timeframe_obs = (obs.T >= start_year & obs.T <= end_year);
        ref_T_years = obs.T(timeframe_obs);
        prcp = obs.var(:,ismember(obs.T, ref_T_years),:); %cru=cru(timeframe_obs);
        prcp_anomaly = prcp - mean(prcp(:,ismember(ref_T_years, anomaly_years),:)); %cru_anomaly = cru - mean(cru);
        prcp_standardized = prcp_anomaly./std(prcp,0,2); %cru_standardized = cru_anomaly/std(cru);
        
	G = load(make_data_filename(variable, start_month, end_month, scenario,'GM'));
        timeframe_m = ismember(single(G.time), ref_T_years);
        anomaly_timeframe = ismember(single(G.time), anomaly_years);
        ref_T_years = G.time(timeframe_m);
        end_year = ref_T_years(end);
        
	GA = load(make_analysis_filename(variable,scenario, start_year, end_year, 500));
        MMM = GA.MMM.MMM(:,timeframe_m,:);                
        MMM_anomaly = MMM - mean(MMM(:,anomaly_timeframe,:),2); 
        MMM_standardized = MMM_anomaly./std(MMM_anomaly,0,2);

        GM = G.GMs(:,timeframe_m,:);
        GM_anomalies = GM - mean(GM(:,anomaly_timeframe,:),2);   
        GM_standardized = GM_anomalies./std(GM_anomalies, 0, 2);

    if(~contains(scenario, {'fast'}))
	I = load(make_data_filename(variable, start_month, end_month, scenario,'all'));
        runs = I.runs(:,timeframe_m,:); 
        runs_anomalies = runs - mean(runs(:,anomaly_timeframe,:),2);   
        runs_standardized = runs_anomalies./std(runs_anomalies, 0, 2);
    end

        figure(1); hold off; clf;

        Indx = size(MMM, 3);
        for dx = 1:Indx
            subplot(2*Indx,2,2*dx-1)
            set(gca,'FontSize',16); %set(gca,'LineWidth',2);
            title('a. Standardized');
            ylabel("Standardized Precipitation Anomaly"); hold on;
            
            plot(ref_T_years,runs_standardized(:,:,dx), '-', 'Color', bckgnd, 'LineWidth', .1, 'HandleVisibility', 'off'); 
            plot(ref_T_years,GM_standardized(:,:,dx), 'Color', mdgnd, 'LineStyle', '-', 'LineWidth', 1, 'HandleVisibility', 'off');
            plot(ref_T_years,MMM_standardized(:,:,dx),'Color', clr, 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off');
            plot(ref_T_years,prcp_standardized(:,:,dx), 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
            %plot(ref_T_years,cru_standardized, '--', 'Color', .8*[1,1,1], 'HandleVisibility', 'off');
            ylim([-4,4]); xlim([start_year, end_year]);

            subplot(2*Indx,2,2*dx); hold on; %yikes again
            set(gca,'FontSize',16); %set(gca,'LineWidth',2);
            title('b. Anomalies');
            ylabel("Precipitation Anomaly (mm/day)"); 
            p_runs_s = plot(ref_T_years,runs_anomalies(:,:,dx), '-', 'Color', bckgnd, 'LineWidth', .1);
            p_gm_s = plot(ref_T_years,GM_anomalies(:,:,dx), '-', 'Color', mdgnd, 'LineWidth', 1);
            p_mmm_s = plot(ref_T_years,MMM_anomaly(:,:,dx),'-', 'Color', clr, 'LineWidth', 2);
            p_actual_s = plot(ref_T_years,prcp_anomaly(:,:,dx), 'k-', 'LineWidth', 2);
            %p_actual_s_cru = plot(ref_T_years,cru_anomaly, '--', 'Color', .8*[1,1,1]);
            xlim([start_year, end_year]); ylim([-2,2]);

            legend([p_runs_s(1), p_gm_s(1), p_mmm_s, p_actual_s], 'Runs', 'IMs', 'MMM', 'GPCC', 'Location', 'south');%, p_actual_s_cru]'northwest'); 'CRU', %
        end
        if(save)
            savefig(['figures/', variable, '/', scenario, '_Fig1p_', num2str(ref_T_years(1)), '-', num2str(ref_T_years(end))]);
        end
    end
end

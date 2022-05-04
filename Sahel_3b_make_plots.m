save = true;
start_year=1901;
anomaly_years = 1901:1950;
%should work with TS if I pick a basin now, but no need to make this figure
%for that.
start_month = 7;%7
end_month = 9;%9

%TODO put colors other than blue for the different scenarios!
%scenarios = {'cmip6_hfast', 'cmip6_afast', 'cmip6_nfast', 'cmip6_gfast'};
scenarios = {'cmip6_h'};%'amip-hist','amip-piF'};

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
        prcp_anomaly = prcp - mean(prcp); %(:,ismember(ref_T_years, anomaly_years) %cru_anomaly = cru - mean(cru);
        prcp_standardized = prcp_anomaly./std(prcp,0,2); %cru_standardized = cru_anomaly/std(cru);
        prcp_smth = smoothdata(prcp_anomaly, 2, 'movmean', 20);
        loc_mean = mean(prcp_smth, 2);
        prcp_smth_s = smoothdata((prcp_smth-loc_mean)./std(prcp_smth, 0,2)+loc_mean,2,'movmean',5);
        prcp_smth = smoothdata(prcp_smth, 2, 'movmean',5);
        
	G = load(make_data_filename(variable, start_month, end_month, scenario,'GM'));
        timeframe_m = ismember(single(G.time), ref_T_years);
        %trying removing this and comparing each model to its clim instead
        anomaly_timeframe = ismember(single(G.time), anomaly_years);
        ref_T_years = G.time(timeframe_m);
        end_year = ref_T_years(end);
        
	GA = load(make_analysis_filename(variable,scenario, start_year, end_year, 500));
        MMM = GA.MMM.MMM(:,timeframe_m,:);                
        MMM_anomaly = MMM - mean(MMM(:,anomaly_timeframe,:),2); 
        MMM_standardized = MMM_anomaly./std(MMM_anomaly,0,2);
        MMM_smth = smooth(MMM_anomaly, 20)';
        loc_mean = mean(MMM_smth, 2);
        MMM_smth_s = smoothdata((MMM_smth-loc_mean)/std(MMM_smth, 0, 2)+loc_mean,2,'movmean',5);
        MMM_smth = smooth(MMM_smth, 5);

        GM = G.GMs(:,timeframe_m,:);
        GM_anomalies = GM - mean(GM(:,anomaly_timeframe,:),2);   
        GM_standardized = GM_anomalies./std(GM_anomalies, 0, 2);
        GM_smth = smoothdata(GM_anomalies, 2, 'movmean', 20);
        loc_mean = mean(GM_smth, 2);        
        GM_smth_s = smoothdata((GM_smth-loc_mean)./std(GM_smth, 0, 2)+loc_mean, 2, 'movmean', 5);
        GM_smth = smoothdata(GM_smth, 2, 'movmean', 5);

    if(~contains(scenario, {'fast'}))
    %individual models don't have mean subtracted yet
	I = load(make_data_filename(variable, start_month, end_month, scenario,'all'));
    M = load(make_data_filename(variable, start_month, end_month, scenario,'MM'));
        I = struct2table(rmfield(I, 'time')); I = I(ismember(I.model(:,1), M.models(:,1)),:);
        runs = I.runs(:,timeframe_m,:); 
        [~,L] = ismember(I.model(:,2), M.models(:,2));
        runs_anomalies = runs - mean(M.MMs(L, anomaly_timeframe,:),2);
        %runs_anomalies = runs - mean(runs(:,anomaly_timeframe,:),2);   
        runs_standardized = runs_anomalies./std(runs_anomalies, 0, 2);
        runs_smth = smoothdata(runs_anomalies, 2, 'movmean', 20);
        loc_mean = mean(runs_smth, 2);        
        runs_smth_s = smoothdata((runs_smth-loc_mean)./std(runs_smth, 0, 2)+loc_mean, 2, 'movmean', 5);
        runs_smth = smoothdata(runs_smth, 2, 'movmean', 5);
        [L,H] = confidence_interval(runs_smth, 1);
        [L_s, H_s] = confidence_interval(runs_smth_s,1);
        sf = @(X) smoothdata(smoothdata(X, 2, 'movmean', 20), 2, 'movmean', 5);
        L = sf(L); H = sf(H); L_s = sf(L_s); H_s = sf(H_s);
    end

        figure(1); hold off; clf;

        Indx = size(MMM, 3);
        for dx = 1:Indx
            subplot(2*Indx,2,2*dx-1)
            set(gca,'FontSize',16); %set(gca,'LineWidth',2);
            title('a. Standardized');
            ylabel("Standardized Precipitation Anomaly"); hold on;
            
            plot(ref_T_years,runs_smth_s(:,:,dx), '-', 'Color', bckgnd, 'LineWidth', .1, 'HandleVisibility', 'off'); 
            plot(ref_T_years, L_s(:,:,dx), 'k:', 'LineWidth',2)
            plot(ref_T_years, H_s(:,:,dx), 'k:', 'LineWidth',2)
            for gm = 1:size(GM_smth_s,1)
                mod_name = G.models{gm};
                if(strcmp(mod_name, GA.indiv_s.best_models_r{1}))
                    cm = [0.47,0.67,0.19];
                    hv = 'on';
                    lw = 2;
                else
                    cm = mdgnd;
                    hv = 'off';
                    lw = 1;
                end
                plot(ref_T_years,GM_smth_s(gm,:,dx), 'Color', cm, 'LineStyle', '-',...
                'LineWidth', lw, 'HandleVisibility', hv, 'DisplayName', mod_name);
            end
            plot(ref_T_years,MMM_smth_s(:,:,dx),'Color', clr, 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off');
            plot(ref_T_years,prcp_smth_s(:,:,dx), 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
            %plot(ref_T_years,cru_standardized, '--', 'Color', .8*[1,1,1], 'HandleVisibility', 'off');
            %ylim([-4,4]); 
            xlim([start_year, end_year]);

            subplot(2*Indx,2,2*dx); hold on; %yikes again
            set(gca,'FontSize',16); %set(gca,'LineWidth',2);
            title('b. Anomalies');
            ylabel("Precipitation Anomaly (mm/day)"); 
            p_runs_s = plot(ref_T_years,runs_smth(:,:,dx), '-', 'Color', bckgnd, 'LineWidth', .1);
            plot(ref_T_years, L(:,:,dx), 'k:', 'LineWidth',2, 'HandleVisibility', 'off')
            plot(ref_T_years, H(:,:,dx), 'k:', 'LineWidth',2, 'HandleVisibility', 'off')
            for gm = 1:size(GM_smth,1)
                mod_name = G.models{gm};
                if(strcmp(mod_name, GA.indiv_s.best_models_e{1}))
                    cm = [0.47,0.67,0.19];
                    hv = 'on';
                    lw=2;
                else
                    cm = mdgnd;
                    hv = 'off';
                    lw=1;
                end
                p_gm_s = plot(ref_T_years,GM_smth(gm,:,dx), '-', 'Color', cm,...
                    'LineWidth', lw,'HandleVisibility', 'off', 'DisplayName', G.models{gm});
            end
            p_mmm_s = plot(ref_T_years,MMM_smth(:,:,dx),'-', 'Color', clr, 'LineWidth', 2);
            p_actual_s = plot(ref_T_years,prcp_smth(:,:,dx), 'k-', 'LineWidth', 2);
            %p_actual_s_cru = plot(ref_T_years,cru_anomaly, '--', 'Color', .8*[1,1,1]);
            xlim([start_year, end_year]); yl = ylim; yl = max(abs(yl)); ylim([-yl,yl]);

            legend([p_runs_s(1), p_gm_s(1), p_mmm_s, p_actual_s], 'Runs', 'IMs', 'MMM', 'GPCC', 'Location', 'southwest');%, p_actual_s_cru]'northwest'); 'CRU', %
        end
        if(save)
            savefig(['figures/', variable, '/', scenario, '_Fig1p_smth_', num2str(ref_T_years(1)), '-', num2str(ref_T_years(end))]);
        end
    end
end

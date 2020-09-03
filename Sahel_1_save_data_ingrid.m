%{
Use for all data accessed without jupyter and the cloud, rather through ingrid. 
Includes CMIP5 and Vanilla AMIP simulations, as well as certain AMIP+RAD
simulations.

Downloads data for all models and runs used in Herman et al. 2020 
and in its followup paper (in progress), saving them in a single file.

Log files for seach scenario "<Scenario>_log.txt" detail which simulations
were skipped and why.
TODO automate the log and the models.mat creation

%CMIP6 institutions file calculated using the following code, after 
downloading CMIP6 historical simulations using jupyter:
H = load('data/pr/cmip6_h_all.mat');
[~, ia, ~] = unique(H.model(:,2));
institutions = H.model(ia, 1:2);
save('data/institutions_cmip6.mat', 'institutions');
%}

%TODO: update the file naming functions in observations
%TODO: just save all available years?

clear
get_observations = false;

%Create MONTH_NAMES and MONTH_DAYS matrices for use in changing units later
month_names = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
month_days  = 31*ones(1,12); month_days(2)=28; %should it just be uniformly 30-day months? I think not in CESM
month_days(ismember(month_names,{'Sep','Apr','Jun','Nov'}))=30;
%{
Simulations = {'Coupled', 'piC', 'AMIP', 'AMIP+RAD'};
Variables = {'pr'; 'ts', 'pr'; 'ts', 'pr', 'pr'};
SST_indices = {'NA', 'GT'};
MIPs = {'CMIP5'; 'CMIP6', 'CMIP5'; 'CMIP6', 'AMIP5'; 'AMIP6pF', 'AMIP6'; 'ERA20CM'; 'NOAA'};
%}
scenarios = {'VanillaAMIP'};%'PSL_FACTS'};%'historical','historicalAerosol','historicalNat','historicalGHG','piControl'};%'historicalMisc' is volcanoes only.  
shortcuts = {'amip-piF'};%'p'};%'h', 'a', 'n', 'g', 'piC'};%, 'a6'};
Generation = 6;
variable = 'pr';

start_month = [7];
end_month = [9];
start_year = 1901;
end_year = 2014;

global ref_T_years
ref_T_years = (start_year:1:end_year)';

s_lat = 12;
n_lat = 18;

umbrella=load(['data/institutions_cmip', num2str(Generation), '.mat']);

space_range_average = [... also includes some background time stuff...
    '/.pr/',...   
    'time//T/renameGRID/',... rename grid if needed.
    'T/1/monthlyAverage/',... changes units from months to years
    'lon/-20/40/RANGE/',... zonal range
    'lat/', num2str(s_lat), '/', num2str(n_lat), '/RANGE',... meridional range
    '%7B/lat/cosd/%7D/',... box size changes with lat
    '%5Blon/lat%5D/weighted-average/']; % zonal and meridional average

for i = 1:length(start_month)
    month0 = month_names{start_month(i)};
    month1 = month_names{end_month(i)};
    mkdir(['data/', month0, '-', month1])
    
    %TODO this is only true if models don't just use 30-day months....
    days_per_month = mean(month_days(start_month(i):end_month(i)));
    
    %If I'm going to do a 9-year running average CENTERED at 1901 and 1999,
    %then I need to start at 1897 and end at 2003.
    %It seems that the observational data doesn't start before 1901!
    time_range = ['T/%28', month0, '%20', num2str(start_year), '%29%28',...
        month1, '%20', num2str(end_year), '%29RANGE/'];
    seasonal_average = ['T/%28', month0, '-', month1,'%29/seasonalAverage/',...
        'T/12/STEP/dods'];
    url_end_h = [space_range_average, time_range, seasonal_average];
    url_end_piC = [space_range_average, seasonal_average];
    
    %Observations
    if(get_observations)
        %TODO these are old observations file and attribute names
        observations_url = [...
            'http://iridl.ldeo.columbia.edu/SOURCES/',...
            '.WCRP/.GCOS/.GPCC/.FDP/.version7/.1p0/.prcp/',...
            'X/-20/40/RANGE/',...
            'Y/',num2str(s_lat), '/', num2str(n_lat), '/RANGE/',...
            '%7B/Y/cosd/%7D/',...
            '%5BX/Y%5D/weighted-average/', time_range, seasonal_average];

        observations_url_CRU = [...
            'http://kage.ldeo.columbia.edu:81/expert/home/.datasets/.CRU3.25/.cru3p25.nc/.pre/',...
            'X/-20/40/RANGE/',...
            'Y/',num2str(s_lat), '/', num2str(n_lat), '/RANGE',...
            '%7B/Y/cosd/%7D/',...
            '%5BX/Y%5D/weighted-average/', time_range, seasonal_average];

        % TODO fix this to a weighted average and put back in
        %SST_forced_url = 'http://iridl.ldeo.columbia.edu/SOURCES/.IRI/.FD/.NSIPP-1/.History/.ensemble/.monthly/.prcp/Y/cosd/mul/X/-20/40/RANGE/Y/12/18/RANGE/%5BX/Y%5Daverage/T/%28Jul%201950%29/%28Sep%202000%29/RANGE/T/%28Jul-Sep%29/seasonalAverage/T/12/STEP/dods';

        prcp_unstandardized = ncread(observations_url, 'prcp'); 
        ref_T_months = ncread(observations_url,'T');
        ref_T_years = adjust_time(ref_T_months, 1960);
        L = length(ref_T_years);

        %TODO: update make observations file name function
        ObsFile = matfile(make_observations_file_name(month0, month1),'Writable', true);
        ObsFile.prcp(1,1:L) = rainfall_to_days(prcp_unstandardized', days_per_month);
        %To avoid redundancy, we only save the years for the observations, and
        %check that the years for the model outputs would have matched.
        ObsFile.T(1,1:L) = ref_T_years';

        CRU_obs = matfile(['data/', month0, '-',month1,'/CRU_data.mat'], 'Writable', true);
        CRU_obs.prcp(1,1:L) = rainfall_to_days(ncread(observations_url_CRU, 'pre')', days_per_month);
    end
    
    for j = 1:length(scenarios)
        clear model
        scenario = scenarios{j};
        sc = shortcuts{j};
        log = fopen(['data/', variable, '/', sc, '_log.txt'], 'wt');
        switch scenario
            case {'historical','historicalAerosol','historicalNat','historicalGHG','piControl'} %cmip5
                url_setup = 'http://strega.ldeo.columbia.edu:81/expert/CMIP5/.byScenario/.';
            case 'amip' %cmip6, many start at 1950. supplemented with ERA20CM (e) and PSL-FACTS (p)
                url_setup = 'http://carney.ldeo.columbia.edu:81/expert/home/.OTHER/.rebecca/.netcdf/.cmip6/.amip/.';
            case 'PSL_FACTS'
                url_setup = 'http://carney.ldeo.columbia.edu:81/expert/home/.OTHER/.rebecca/.netcdf/PSL-FACTS/.';
            case 'VanillaAMIP'
                url_setup = 'http://carney.ldeo.columbia.edu:81/expert/home/.OTHER/.rebecca/.netcdf/.cmip6/.amip-piForcing/.';
        end

        model_file_name = ['data/',variable, '/', sc,'_all.mat'];
        if(exist(model_file_name,'file')==2)
            load(model_file_name);
        end
        fprintf("Accessing scenario %s\n", scenario);
        fprintf("Opening file %s_models.txt\n", scenario);
        models = fopen([scenario, '_models.txt']);
        
        file_name = fgetl(models);
        next_line=1;
        while ischar(file_name)
            elements = strsplit(file_name, '_');
            tokeep = false;
            model_name = get_model_name(elements, scenario);
            [mpname, run_name] = make_model_and_p_name(elements, scenario); %already type char
            if(exist('model','var') && any(contains(model(:,2), mpname) & contains(model(:,3), run_name)))
                next_line=next_line+1;
                fprintf(log, "Already saved run %s\n", file_name);
                continue
            end
            if(strcmp(scenario, "piControl"))
                url = [url_setup, file_name, url_end_piC];
            else
                url = [url_setup, file_name, url_end_h];
            end
            T_units = ncreadatt(url, 'T', 'units');
            P_units = ncreadatt(url, 'pr', 'units');
            if(check_units(T_units, P_units, log))
                non_standardized_run = ncread(url, 'pr');
                l = length(non_standardized_run);
                T=ncread(url, 'T');
                T_years = adjust_time(T, get_ref_year(T_units));
                %Why did I comment this out? Because the amip runs
                %were short?
                if(~check_years(scenario, T_years, l, model_name, run_name, log))
                    continue
                end
                if(sum(isnan(non_standardized_run))>0)
                    fprintf(log, "%s contains NaN in run %s\n",...
                        model_name, run_name);
                else
                    model(next_line, 2:3) = {mpname, run_name};
                    runs(next_line,1:l)   = flux_to_rainfall(non_standardized_run');
                    time(next_line,1:l)   = T_years;
                    next_line=next_line+1;
                    tokeep = true;
                    save(model_file_name,'model','runs', 'time');
                end
            end
            if(tokeep)
                fprintf(log, "Saved run %s with %u years.\n", file_name, l);
            else
                fprintf(log, "skipping file %s\n", file_name);
            end
            file_name = fgetl(models);
        end
        fclose(models);
        full_model_names = model(:,2); %just the full model names
        umbrella_model_names = find_umbrella_names(full_model_names, umbrella);
        model(1:next_line-1,1) = umbrella_model_names;
        save(model_file_name, 'model', 'runs', 'time');
        clear model runs
    end
end

%functions for legibility
function units_ok = check_units(T_units, P_units, log)
    units_ok = all([(T_units(1:13) == 'months since '),...
        (T_units(18:end)=='-01-01'),...
        (P_units=='kg m-2 s-1')]);
    if(~units_ok)
        fprintf(log, "cannot comprehend units: %s\n",T_units);
    end
end

function years_ok = check_years(scenario, T_years, l, model_name, run, log)
    global ref_T_years;
    L = length(ref_T_years);
    %We don't want to check the years for the pre-industrial
    %controls, because the years are meaningless. We keep ALL
    %AVAILABLE YEARS -- thus the length might be more than 107!
    %It seems that we can gather this data even though we put
    %years into the URL (but I should double-check!) since the
    %years in the URL don't overlap with the available years.
    %But what if some are less than 107 but equal to 100?
    %perhaps I still don't want to throw them out?
    if(strcmp(scenario, "piControl"))
        years_ok = (L<=l);
        if(~years_ok)
            fprintf(log, "incomplete data for model %s, run %s: %u years\n",...
                model_name, run, l);
        end
    else
        %the fact that I used a floor here is not
        %consistent with my subtracting .5 from the years
        %before. Actually maybe it's not, because the
        %offset will always be more than half a year...
        years_ok = (length(T_years) == L && all(T_years == ref_T_years));
        if(~years_ok) 
            fprintf(log, "incomplete data for model %s, run %s: %u:%u\n",...
                model_name, run, T_years(1),T_years(length(T_years)));
        end
    end
end

function ref_year = get_ref_year(T_units)
    ref_year = str2num(T_units(14:17));
end

%IS THIS EVEN RIGHT??? OR IS IT A WEIGHTED MEAN???
%mm/month*month/avg#days=mm/day
function [flux] = rainfall_to_days(rainfall, avg_days_per_month)
    flux = rainfall/(avg_days_per_month);
end

%data are reported in mid-months. We take the floor so that we get the data
%point centered on the year.
function [T] = adjust_time(Tm, refT)
    T = floor(Tm/12+refT);
end

%Extracts name from line with name and runs
function name = get_model_name(line, scenario)
    % special version for PSL-facts
    if(strcmp(scenario, 'PSL_FACTS'))
        name = line{:,2};
    else
        name = line{:,3};
    end
end

%Concatenates name with p number so that distinct p numbers will be saved as separate files.
%If we want to collapse all offshoots of one model into one, we would want 
%to average over these. But that would happen in an analysis file.
%TODO: changes here to include f can be general. dif between split and
%strsplit seems to just be the dimension of the output array.
function [fname, run] = make_model_and_p_name(line, scenario)
    if(strcmp(scenario, 'PSL_FACTS'))
        run = line{6};
        run_stats = split(run, ["s","."]);
        fname = get_model_name(line, scenario);
        run = run_stats{2};
    else
        run = line{5};    
        run_stats = split(run, 'p');
        p_and_f = split(run_stats(2), 'f');
        p = p_and_f{1}; 
        fname = [get_model_name(line, scenario), ' p', p];
    end
end

function[umbrella_names] = find_umbrella_names(used_models, umbrella)
   umbrella_names = cell(length(used_models),1);
   %case for cmip6
   if(~isfield(umbrella, 'abbrev'))
       [~,Loc] = ismember(used_models,umbrella.institutions(:,2));
       umbrella_names = umbrella.institutions(Loc,1);
   %case for cmip5
   else
       for k = 1:length(umbrella.models)
           umbrella_names(contains(used_models, umbrella.models(k)), 1) = umbrella.abbrev(k);
       end
   end
end

%mm/month*month/avg#days*day/24hours*hour/60min*min/60s*m/1000mm*1000kg/m^3=kg/s/m^2
%kg/s/m^2*m^3/1000kg*1000mm/m*60s/min*60min/hour*24hours/day = mm/day
function [flux] = flux_to_rainfall(rainfall)
    flux = rainfall*24*60*60;
end

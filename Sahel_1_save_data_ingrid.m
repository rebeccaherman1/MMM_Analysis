%{
Use for all data accessed without jupyter and the cloud, rather through the IRI library.
Data is preprocessed using ingrid (see http://rainbow.ldeo.columbia.edu/documentation/libraries/ingrid/Ingrid_5.html#SEC39). 
Includes CMIP5 and observations. Can also be used for short CMIP5 AMIP simulations.

Log files for seach scenario "<Scenario>_log.txt" detail which simulations
were skipped and why.
%}

clear

%customizable variables
get_observations = true;
get_simulations = false;
scenarios = {'historical','historicalAerosol','historicalNat','historicalGHG','piControl'};%'historicalMisc' is volcanoes only.  
%there are additional amip simulations which I downloaded onto a Lamont
%server which can be accessed here, called 'amip', 'PSL-FACTS', and
%'VanillaAMIP'. They are not included here because they were too short and
%redundant to CMIP6 amip-piF and amip-hist runs. Be sure to use the cmip6 
%institution file for those.
shortcuts = {'h', 'a', 'n', 'g', 'piC'};
%make sure the shortcuts match the scenarios if they are modified.
variable = 'ts';
start_month = 7;
end_month = 9;

%define area limits
switch variable
    case 'pr'
        s_lat = 12;
        n_lat = 18;
        e_lon = -20;
        w_lon = 40;
    case 'ts'
        %NA
        s_lat = 10;
        n_lat = 40; 
        e_lon = -75;
        w_lon = -15; 
        %GT
        GT_lat = [-20,20]; 
end

%create data folder
drnm = ['data/', variable];
if(~exist(drnm, 'dir'))
    mkdir(drnm)
end

%Create MONTH_NAMES and MONTH_DAYS matrices for use in changing units later
month_names = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
month_days  = 31*ones(1,12); month_days(2)=28; %should it just be uniformly 30-day months? I think not in CESM
month_days(ismember(month_names,{'Sep','Apr','Jun','Nov'}))=30;

%ingrid preparations
month0 = month_names{start_month};
month1 = month_names{end_month};
%NOTE: this is only true if models don't just use 30-day months....
days_per_month = mean(month_days(start_month:end_month));
seasonal_average = ['T/%28', month0, '-', month1,'%29/seasonalAverage/',...
    'T/12/STEP/dods'];
space_range_average = [... 
    %'/time//T/renameGRID',... rename grid if needed. causes cold to fail
    %if not needed, hence it is commented out.
    '/T/1/monthlyAverage',... changes units from months to years
    '/lon/', num2str(e_lon), '/', num2str(w_lon), '/RANGE',... zonal range
    '/lat/', num2str(s_lat), '/', num2str(n_lat), '/RANGE',... meridional range
    '%7B/lat/cosd/%7D/',... box size changes with lat
    '%5Blon/lat%5D/weighted-average/']; % zonal and meridional average
url_end = [space_range_average, seasonal_average]; 
if(strcmp(variable, 'ts'))
    GT_space_average = [... 
                        '/time//T/renameGRID',... rename grid if needed.
                        '/T/1/monthlyAverage',... changes units from months to years
                        '/lat/', num2str(GT_lat(1)), '/', num2str(GT_lat(2)), '/RANGE',... meridional range
                        '%7B/lat/cosd/%7D/',... box size changes with lat
                        '%5Blon/lat%5D/weighted-average/']; % zonal and meridional average
    land_mask = [... %used for GT in simulations where ts includes land measurements
                'SOURCES/.WORLDBATH/.bath/',... 
                'X/%28lon%29renameGRID/',...
                'Y/%28lat%29renameGRID',...
                '%5Blon/lat%5DregridLinear/',...
                '0/maskgt/0/mul/add/'];
    url_end_GT = [land_mask, GT_space_average, seasonal_average];
end

%Observations
if(get_observations)
    obs_file_name = make_data_filename(variable, start_month, end_month, 'observations', 0);
    fprintf('Writing file %s\n', obs_file_name);
    ObsFile = matfile(obs_file_name,'Writable', true);
    %matfile(['data/', variable, '/', num2str(start_month), '-', num2str(end_month), '/observations.mat']
    switch variable
        case 'pr'
            source = '.WCRP/.GCOS/.GPCC/.FDP/.version2018/.1p0/.prcp/';
            ObsFile.src = 'GPCC version18';
            v = 'prcp';
        case 'ts'
            source = '.NOAA/.NCDC/.ERSST/.version5/.sst/';
            ObsFile.src = 'ERSSTv5';
            v = 'sst';
    end
    obs_url_setup = [...
        'http://iridl.ldeo.columbia.edu/SOURCES/',...
        source];
    observations_url_space = [...
        'X/', num2str(e_lon), '/', num2str(w_lon), '/RANGE/',...
        'Y/',num2str(s_lat), '/', num2str(n_lat), '/RANGE/'];
    obs_url_averaging = [...
        '%7B/Y/cosd/%7D/',...
        '%5BX/Y%5D/weighted-average/', seasonal_average];
    observations_url = [obs_url_setup, observations_url_space,obs_url_averaging];
    obs_units = ncreadatt(observations_url, 'T', 'units');
    var_unstandardized = ncread(observations_url, v);
    if(strcmp(variable, 'pr'))
        var_unstandardized = rainfall_to_days(var_unstandardized, days_per_month);
    end
    ObsFile.var = var_unstandardized;
    ObsFile.T = adjust_time(ncread(observations_url,'T'), get_ref_year(obs_units))';
    
    if strcmp(variable, 'pr') %add an additional file for CRU observations
        observations_url_CRU = [...
            'http://kage.ldeo.columbia.edu:81/expert/home/.datasets/.CRU3.25/.cru3p25.nc/.pre/',...
            observations_url_space,obs_url_averaging];
        
        CRU_obs = matfile(['data/', variable,'/CRU_data.mat'], 'Writable', true);
        CRU_obs.var = rainfall_to_days(ncread(observations_url_CRU, 'pre')', days_per_month);
        CRU_obs.T = adjust_time(ncread(observations_url_CRU, 'T'), get_ref_year(ncreadatt(observations_url_CRU, 'T', 'units')));
        CRU_obs.src = 'CRU';
    else %add additional ocean basins to the observaitons file
        obs_url_GT = [obs_url_setup,... %don't need a land mask for SST observations
            'Y/',num2str(GT_lat(1)), '/', num2str(GT_lat(2)), '/RANGE/',...
            obs_url_averaging];
        var_unstandardized = [var_unstandardized; ncread(obs_url_GT, v)];
        var_unstandardized(3,:) = var_unstandardized(1,:) - var_unstandardized(2,:);
        ObsFile.var = permute(var_unstandardized, [3,2,1]);
        ObsFile.indices = permute({'NA'; 'GT'; 'NARI'},[3,2,1]);
    end
end

%Historical Simulations
if(get_simulations)
    umbrella=load('data/institutions_cmip5.mat');
    for j = 1:length(scenarios)
        clear model
        scenario = scenarios{j};
        sc = shortcuts{j};
        log = fopen(['data/', variable, '/', sc, '_log.txt'], 'wt');
        switch scenario
            case {'historical','historicalAerosol','historicalNat','historicalGHG','piControl'} %cmip5
                url_setup = ['http://strega.ldeo.columbia.edu:81/expert/CMIP5/.byScenario/.', scenario,'/.atmos/.mon/.', variable, '/.'];
            case 'amip' %cmip6, many start at 1950. supplemented with ERA20CM (e) and PSL-FACTS (p)
                url_setup = 'http://carney.ldeo.columbia.edu:81/expert/home/.OTHER/.rebecca/.netcdf/.cmip6/.amip/.';
            case 'PSL_FACTS'
                url_setup = 'http://carney.ldeo.columbia.edu:81/expert/home/.OTHER/.rebecca/.netcdf/PSL-FACTS/.';
            case 'VanillaAMIP'
                url_setup = 'http://carney.ldeo.columbia.edu:81/expert/home/.OTHER/.rebecca/.netcdf/.cmip6/.amip-piForcing/.';
        end

	model_file_name = make_data_filename(variable, start_month, end_month, sc,'all');
        if(exist(model_file_name,'file')==2)
            load(model_file_name);
        end
        fprintf("Accessing scenario %s\n", scenario);
        fprintf("Opening file %s_models.txt\n", scenario);
        models = fopen(['data/',variable, '/',scenario, '_models.txt']);

        file_name = fgetl(models);
        next_line=1;
        %for CMIP5, one line corresponds to all runs from one model. For
        %data downloaded directly from ESGF, one line corresponds to one
        %run.
        while ischar(file_name)
            elements = strsplit(file_name, {'_', '\t', ' '});
            model_name = get_model_name(elements, scenario);
            runs_ids = get_runs(elements, scenario);
            for rx = 1:length(runs_ids)
                r = runs_ids{rx};
                tokeep = false;
                [mpname, run_name] = make_model_and_p_name({model_name, r}, scenario); %already type char
                if(exist('model','var') && any(contains(model(:,2), mpname) & contains(model(:,3), run_name)))
                    next_line=next_line+1;
                    fprintf(log, "Already saved run %s\n", file_name);
                    continue
                end
                url = [url_setup,model_name, '/.',r, '/.', variable, url_end];
                T_units = ncreadatt(url, 'T', 'units');
                P_units = ncreadatt(url, variable, 'units');
                if(check_units(T_units, P_units, variable, log))
                    non_standardized_run = ncread(url, variable);
                    l = length(non_standardized_run);
                    T=ncread(url, 'T');
                    T_years = adjust_time(T, get_ref_year(T_units));
                    if(~check_years(scenario, T_years, l, model_name, run_name, log))
                        continue
                    end
                    if(sum(isnan(non_standardized_run))>0)
                        fprintf(log, "%s contains NaN in run %s\n",...
                            model_name, run_name);
                    else
                        model(next_line, 2:3) = {mpname, run_name};
                        if strcmp(variable, 'pr')
                            non_standardized_run = flux_to_rainfall(non_standardized_run)';
                            s3=1;
                        else
                            url_GT = [url_setup,model_name, '/.',r, '/.', variable, url_end_GT];
                            non_standardized_run = [non_standardized_run, ncread(url, variable)];
                            non_standardized_run(:,3) = non_standardized_run(:,1)-non_standardized_run(:,2);
                            non_standardized_run = permute(non_standardized_run, [3,1,2]);
                            s3=3;
                        end
                        runs(next_line,1:l,1:s3)   = non_standardized_run;
                        time(next_line,1:l)   = T_years;
                        next_line=next_line+1;
                        tokeep = true;
                        save(model_file_name,'model','runs', 'time');
                    end
                end
                if(tokeep)
                    fprintf(log, "Saved model %s run %s with %u years.\n", model_name, run_name, l);
                else
                    fprintf(log, "skipping file %s\n", file_name);
                end
            end
            file_name = fgetl(models);
        end
        fclose(models);
        umbrella_model_names = find_umbrella_names(model_name, umbrella);
        model(1:next_line-1,1) = umbrella_model_names;
        save(model_file_name, 'model', 'runs', 'time');
        clear model runs
    end
    if(strcmp(variable, 'ts'))
        load(model_file_name);
        indices = permute({'NA'; 'GT'; 'NARI'}, [3,2,1]);
        save(model_file_name, 'model', 'runs', 'time', 'indices');
    end
end

%functions for legibility
function units_ok = check_units(T_units, P_units, variable, log)
    T_units_ok = all([(T_units(1:13) == 'months since '),...
                    (T_units(18:end)=='-01-01')]);
    var_units_ok = strcmp(variable, 'pr') && strcmp(P_units,'kg m-2 s-1') || ...
        strcmp(variable, 'ts') && strcmp(P_units, 'Kelvin_scale');
    units_ok = T_units_ok && var_units_ok;
    if(~units_ok)
        fprintf(log, "cannot comprehend units: %s\n",T_units);
    end
end

function years_ok = check_years(scenario, T_years, l, model_name, run, log)
    ref_T_years = 1901:2003; %year limits for CMIP5.
    L = length(ref_T_years);
    %For piC, just make sure there are enough years total.
    if(strcmp(scenario, "piControl"))
        years_ok = (L<=l);
        if(~years_ok)
            fprintf(log, "incomplete data for model %s, run %s: %u years\n",...
                model_name, run, l);
        end
    else
        %for historical simulations, check that the years include all the necessary ones. 
        years_ok = sum(T_years >= min(ref_T_years) & T_years <= max(ref_T_years)) == L;
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
    %files downloaded directly from ESGF are listed by name individually
    if(strcmp(scenario, 'PSL_FACTS') || strcmp(scenario, 'VanillaAMIP'))
        name = line{2};
    else %for CMIP5, the runs are listed within one line after the model name.
        name = line{1};
    end
end

function runs = get_runs(line, scenario)
    if(strcmp(scenario, 'PSL_FACTS') || strcmp(scenario, 'VanillaAMIP'))
        runs = line(6);
    else
        runs = line(2:end);
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
        run = line{2};    
        run_stats = split(run, 'p');
        %useful only for post-CMIP5 simulations
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

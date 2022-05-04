%saves weighted MMs (model means) and IMs (institution means) as 
%"<scenario>_XM.mat". Weights to use for the following tier are stored in 
%TRUST. 
%uses only models which contribute AA simulations.

%TODO: if I remove the year constraints in 1_save_data, my year checking
%here might not work...

tosave = true;
start_month = 7;
end_month = 9;
gts = false;

realm = 'cmip5';
switch realm
    case 'amip'
        scenarios = {'amip-piF', 'amip-hist'};
        variables = {'pr'};
    case 'cmip5'
        scenarios = {'h', 'a', 'n', 'g'};
        piCs = 'piC';
        variables = {'pr','ts'}; %, 
    case 'cmip6'
        scenarios = {'cmip6_h','cmip6_a', 'cmip6_n', 'cmip6_g'};
        piCs = 'cmip6_piC';
        variables = {'pr','ts'}; 
    otherwise
        fprintf("what do you want?")
end
CM = scenarios{2};
if(gts)
    variables = {'globalts'};
    scenarios = scenarios(1);
end   

%TODO decide if I want nanmean or not
vert_mean = @(X) mean(X,1); vert_sum = @(X) sum(X,1);
vert_nan_mean = @(X) nanmean(X,1); vert_nan_sum = @(X) nansum(X,1);

for v = 1:length(variables)
    var = variables{v};
    if(strcmp(var, 'globalts'))
    	AA = load(make_data_filename('pr', 7,9, CM, 'all'));%'MM'));%
    else
        AA = load(make_data_filename(var, start_month, end_month, CM, 'all'));%'a', 'MM'));%
    end
    %AA = load(make_data_filename('ts', start_month, end_month, 'cmip6_h', 'all'));%'a', 'MM'));%
    common_models = unique(AA.model(:,1)); %just model for CMIP6?
    for j = 1:length(scenarios)
        clear MM GM
        scenario = scenarios{j};
        fprintf("Accessing scenario %s variable %s\n", scenario, var);

        fn = make_data_filename(var, start_month, end_month, scenario, 'all');
        h = load(fn);
        if(isfield(h, 'indices'))
            h_indices = h.indices;
        end
        if(isfield(h, 'lat'))
            lat = h.lat;
            lon = h.lon;
        end
        h = table(h.model, h.runs, repmat(h.time, size(h.model, 1),1), 'VariableNames', {'model', 'runs', 'time'});
        h = h(ismember(h.model(:,1), common_models),:);
        [model_names, I, model_groupings] = unique(h.model(:,2)); nMM = max(model_groupings);
        num_runs = histcounts(model_groupings(~all(all(isnan(h.runs),2),3)), (0:nMM)+.5)';
        MM.MMs = splitapply(vert_nan_mean, h.runs, model_groupings);
        MM.models = [h.model(I,1), model_names];
        MM.trust = sqrt(num_runs);
        if(any(h.time~=h.time(1,:)))
            fprintf('time does not line up!')
        end

        if(tosave) 
            fname = make_data_filename(var, start_month, end_month, scenario, 'MM');
            delete(fname);
            File = matfile(fname, 'Writable', true);
            fprintf("Writing file %s\n", fname);
            File.MMs = MM.MMs; 
            File.models = MM.models; 
            File.trust = MM.trust; 
            NanSims = h.model(any(any(isnan(h.runs),2),3),:);
            File.NanSims = NanSims;
            File.time = h.time(1,:);
            if(exist('h_indices', 'var'))
                File.indices = h_indices;
            end
            if(isfield('h', 'lat'))
                File.lat = h.lat;
            File.lon = h.lon;
            end
        end
	
        if(~strcmp(realm, 'amip') && ~contains(var, 'global'))
            piC = load(make_data_filename(var, start_month, end_month, piCs, 'all')); 
            piC.runs(piC.runs==0)=NaN; piC_lengths = sum(~isnan(piC.runs(:,:,1)), 2);
            
            if(any(strcmp(piC.model(:,1), 'CESM')))
                piC.model(strcmp(piC.model(:,1), 'CESM'),1)={'NCAR'};
            end
            if(any(strcmp(piC.model(:,1), 'CanESM')))
                piC.model(strcmp(piC.model(:,1), 'CanESM'),1)={'CCCma'};
            end
            if(any(strcmp(piC.model(:,1), 'GISS')))
                piC.model(strcmp(piC.model(:,1), 'GISS'),1)={'NASA'};
            end
            if(any(strcmp(piC.model(:,1), 'NorESM')))
                piC.model(strcmp(piC.model(:,1), 'NorESM'),1)={'Nor'};
            end
            
            %T_piC = table(piC.model, piC.runs, piC.time, piC_lengths, 'VariableNames', {'model', 'runs', 'time', 'length'});
            T_piC = table(piC.model, piC.runs, piC_lengths, 'VariableNames', {'model', 'runs', 'length'});
            if(any(strcmp(T_piC.model(:,1), 'CESM')))
                T_piC.model(strcmp(T_piC.model(:,1), 'CESM'),1)={'NCAR'};
            end
            if(any(strcmp(T_piC.model(:,1), 'CanESM')))
                T_piC.model(strcmp(T_piC.model(:,1), 'CanESM'),1)={'CCCma'};
            end
            if(any(strcmp(T_piC.model(:,1), 'GISS')))
                T_piC.model(strcmp(T_piC.model(:,1), 'GISS'),1)={'NASA'};
            end
            if(any(strcmp(T_piC.model(:,1), 'NorESM')))
                T_piC.model(strcmp(T_piC.model(:,1), 'NorESM'),1)={'Nor'};
            end
            %CHANGED THIS
            if(strcmp(realm, 'cmip5'))
                npcm = 1; M_N = MM.models;
            else
                npcm = 2; M_N = model_names;
            end
            relevant_pC_models = ismember(T_piC.model(:,npcm),M_N);
            T_piC = T_piC(relevant_pC_models, :); 
            %before I do this! let me pretend I have more runs!!!!!
            [Lia, Locb] = ismember(T_piC.model(:,npcm), model_names);
            Locb2 = Locb; Locb2(~Lia)=1; %momentarily lose the 0s so that matlab doesn't crash
            num_pC_reps = num_runs(Locb2); num_pC_reps(~Lia)=1; %replace those numbers with just one repetition
            T_piC = repelem(T_piC, num_pC_reps, 1);
            offsets = arrayfun(@randi, T_piC.length);
            %now I will scramble them a bit!
            for row = 1:size(T_piC.model, 1)
                lr = T_piC.length(row); or = offsets(row);
                if(or>1)
                    T_piC.runs(row, 1:lr) = [T_piC.runs(row, or:lr), T_piC.runs(row, 1:or-1)];
                end
                hts = size(h.time,2);
                if(lr<hts)
                    fprintf('Extending piC simulation %s %s\n',T_piC.model{row,2}, T_piC.model{row,3})
                    T_piC.runs(row, (lr+1):hts) = T_piC.runs(row, 1:(hts-lr));
                end
            end            
            [model_names, I, model_groupings] = unique(T_piC.model(:,2)); nMM = max(model_groupings);
            MM.piC_MMs = splitapply(vert_mean, T_piC.runs, model_groupings);
            MM.piC_models = [T_piC.model(I,1), model_names]; 
            MM.piC_trust = sqrt(histcounts(model_groupings, (0:nMM)+.5)');
            
            %remove clim values. No need to take anomalies of sims ever
            %again!
            [has_piC, usd] = ismember(MM.models(:,2), MM.piC_models(:,2));
            MM.clim(has_piC,:,:) = nanmean(MM.piC_MMs(usd(has_piC), :), 2);
            %for models with no piC simulation, estimate clim using an
            %average.
            if(any(~has_piC))
                leftovers = MM.models(~has_piC,1);
                for l = 1:length(leftovers)
                    leftovers{l} = nanmean(MM.piC_MMs(strcmp(MM.piC_models(:,1), leftovers{l}),:),[1,2]);
                end
                MM.clim(~has_piC,:,:) = cell2mat(leftovers);
            end
            MM.MMs = MM.MMs - MM.clim;
            
            if(tosave)  
                File.piC_MMs= MM.piC_MMs;
                File.piC_models = MM.piC_models; 
                File.piC_trust = MM.piC_trust;
                File.clim = MM.clim;
            end
        end
	
        %skipping MM
        %fname = make_data_filename(var, start_month, end_month, scenario, 'MM');
	%MM = load(fname);

        [GM.models, ~, model_groupings] = unique(MM.models(:,1)); nGM = max(model_groupings);
        weights=splitapply(vert_sum, MM.trust, model_groupings);
        GM.GMs = splitapply(vert_sum, MM.trust.*MM.MMs./weights(model_groupings), model_groupings);
        GM.trust = splitapply(@sum, MM.trust, model_groupings)./sqrt(histcounts(model_groupings, (0:nGM)+.5)');
        GM.MMM = sum(GM.trust.*GM.GMs, 1)/sum(GM.trust);
        if(tosave)
            fname = make_data_filename(var, start_month, end_month, scenario, 'GM');
            delete(fname);
            File = matfile(fname, 'Writable', true);
            fprintf("Writing file %s\n", fname);
            File.GMs= GM.GMs; 
            File.models = GM.models;
            File.trust  = GM.trust; 
            File.MMM = GM.MMM;
            File.time=h.time(1,:);
            if(exist('h_indices', 'var'))
                File.indices = h_indices;
            end
            if(exist('lat', 'var'))
                File.lat = lat;
                File.lon = lon;
            end
        end
	
        if(~strcmp(realm, 'amip') && ~any(contains(variables, 'global')))
            [GM.piC_models, I, model_groupings] = unique(MM.piC_models(:,1)); nGM = max(model_groupings);
            weights=splitapply(vert_sum, MM.piC_trust, model_groupings);
            GM.piC_GMs = splitapply(vert_sum, MM.piC_trust.*MM.piC_MMs./weights(model_groupings), model_groupings);
            GM.piC_trust = splitapply(@sum, MM.piC_trust, model_groupings)./sqrt(histcounts(model_groupings, (0:nGM)+.5)');
            if(tosave)
                File.piC_GMs= GM.piC_GMs;
                File.piC_models= GM.piC_models;
                File.piC_trust = GM.piC_trust;
            end
        end
    end
end

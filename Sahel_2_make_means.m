%saves weighted MMs (model means) and IMs (institution means) as 
%"<scenario>_XM.mat". Weights to use for the following tier are stored in 
%TRUST. 
%uses only models which contribute AA simulations.

%TODO: if I remove the year constraints in 1_save_data, my year checking
%here might not work...

tosave = true;
start_month = 7;
end_month = 9;

realm = 'cmip6';
switch realm
    case 'amip'
        scenarios = {'amip-piF', 'amip-hist'};
        variables = {'pr'};
    case 'cmip5'
        scenarios = {'h', 'a', 'n', 'g'};
        piCs = 'piC';
        variables = {'pr', 'ts'}; 
    case 'cmip6'
        scenarios = {'cmip6_h'};%,'cmip6_a', 'cmip6_n', 'cmip6_g'};
        piCs = 'cmip6_piC';
        variables = {'globalts'};%'pr', 
    otherwise
        fprintf("what do you want?")
end
vert_mean = @(X) mean(X,1); vert_sum = @(X) sum(X,1);

for v = 1:length(variables)
    var = variables{v};
    %AA = load(make_data_filename(var, start_month, end_month, scenarios{2}, 'all'));%'a', 'MM'));%
    AA = load(make_data_filename('ts', start_month, end_month, 'cmip6_h', 'all'));%'a', 'MM'));%
    common_models = unique(AA.model(:,1)); %just model for CMIP6?
    for j = 1:length(scenarios)
        scenario = scenarios{j};
        fprintf("Accessing scenario %s variable %s\n", scenario, var);
	    h = load(make_data_filename(var, start_month, end_month, scenario, 'all'));
        if(isfield(h, 'indices'))
            h_indices = h.indices;
        end
        h = table(h.model, h.runs, repmat(h.time, size(h.model, 1)), 'VariableNames', {'model', 'runs', 'time'});
        h = h(ismember(h.model(:,1), common_models),:);
        [model_names, I, model_groupings] = unique(h.model(:,2)); nMM = max(model_groupings);
        num_runs = histcounts(model_groupings(~any(any(isnan(h.runs),2),3)), (0:nMM)+.5)';
        MM.MMs = splitapply(vert_mean, h.runs, model_groupings);
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
            if(exist('indices', 'var'))
                File.indices = h_indices;
            end
	    if(isfield('h', 'lat'))
	        File.lat = h.lat;
		File.lon = h.lon;
	    end
        end

        if(~strcmp(realm, 'amip') & ~contains(var, 'global'))
	    piC = load(make_data_filename(var, start_month, end_month, piCs, 'all')); piC.runs(piC.runs==0)=NaN; piC_lengths = sum(~isnan(piC.runs(:,:,1)), 2);
            %T_piC = table(piC.model, piC.runs, piC.time, piC_lengths, 'VariableNames', {'model', 'runs', 'time', 'length'});
            T_piC = table(piC.model, piC.runs, piC_lengths, 'VariableNames', {'model', 'runs', 'length'});
            relevant_pC_models = ismember(T_piC.model(:,1),h.model(:,1));
            T_piC = T_piC(relevant_pC_models, :); 
            %before I do this! let me pretend I have more runs!!!!!
            [Lia, Locb] = ismember(T_piC.model(:,2), model_names);
            Locb2 = Locb; Locb2(~Lia)=1; %momentarily lose the 0s so that matlab doesn't crash
            num_pC_reps = num_runs(Locb2); num_pC_reps(~Lia)=1; %replace those numbers with just one repetition
            T_piC = repelem(T_piC, num_pC_reps, 1);
            offsets = arrayfun(@randi, T_piC.length);
            %now I will scramble them a bit!
            for row = 1:size(T_piC.model, 1)
                if(offsets(row)>1)
                    T_piC.runs(row, 1:T_piC.length(row)) = [T_piC.runs(row, offsets(row):T_piC.length(row)), T_piC.runs(row, 1:offsets(row)-1)];
                end
            end            
            [model_names, I, model_groupings] = unique(T_piC.model(:,2)); nMM = max(model_groupings);
            MM.piC_MMs = splitapply(vert_mean, T_piC.runs, model_groupings);
            MM.piC_models = [T_piC.model(I,1), model_names]; 
            MM.piC_trust = sqrt(histcounts(model_groupings, (0:nMM)+.5)');
            if(tosave)  
                File.piC_MMs= MM.piC_MMs;
                File.piC_models = MM.piC_models; 
                File.piC_trust = MM.piC_trust;
            end
        end
        [GM.models, ~, model_groupings] = unique(MM.models(:,1)); nGM = max(model_groupings);
        weights=splitapply(vert_sum, MM.trust, model_groupings);
        GM.GMs = splitapply(vert_sum, MM.trust.*MM.MMs./weights(model_groupings), model_groupings);
        GM.trust = splitapply(@sum, MM.trust, model_groupings)./sqrt(histcounts(model_groupings, (0:nGM)+.5)');
        if(tosave)
	    fname = make_data_filename(var, start_month, end_month, scenario, 'GM');
            delete(fname);
            File = matfile(fname, 'Writable', true);
            fprintf("Writing file %s\n", fname);
            File.GMs= GM.GMs; 
            File.models = GM.models;
            File.trust  = GM.trust; 
            %File.MMM = GM.MMM;
            File.time=h.time(1,:);
            if(exist('h_indices', 'var'))
                File.indices = h_indices;
            end
	    if(isfield('h', 'lat'))
	        File.lat = h.lat;
		File.lon = h.lon;
	    end
        end
        if(~strcmp(realm, 'amip') & ~contains(Var, 'global'))
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

%saves weighted MMs (model means) and IMs (institution means) as 
%"<scenario>_XM.mat". Weights to use for the following tier are stored in 
%TRUST.

%TODO: if I remove the year constraints in 1_save_data, my year checking
%here might not work...
tosave = true;

realm = 'cmip5';
switch realm
    case 'amip'
        scenarios = {'amip-piF', 'amip-hist'};
        variables = {'pr'};
    case 'cmip5'
        scenarios = {'h', 'a', 'n', 'g'};
        piCs = 'piC';
        variables = {'pr', 'ts'}; 
    case 'cmip6'
        scenarios = {'cmip6_h','cmip6_a', 'cmip6_n', 'cmip6_g'};
        piCs = 'cmip6_piC';
        variables = {'pr', 'ts'};
    otherwise
        fprintf("what do you want?")
end
vert_mean = @(X) mean(X,1); vert_sum = @(X) sum(X,1);

for v = 1:length(variables)
    var = variables{v};
    for j = 1:length(scenarios)
        scenario = scenarios{j};
        fprintf("Accessing scenario %s variable %s\n", scenario, var);
        h = load(['data/', var, '/', scenario, '_all.mat']); 
        [~,s2,s3] = size(h.runs);
        [model_names, I, model_groupings] = unique(h.model(:,2)); nMM = max(model_groupings);
        num_runs = histcounts(model_groupings(~any(any(isnan(h.runs),2),3)), 1:nMM);
        MM.MMs = splitapply(vert_mean, h.runs, model_groupings);
        MM.models = [h.model(I,1), model_names];
        MM.trust = sqrt(num_runs);
        if(any(h.time~=h.time(1,:)))
            fprintf('time does not line up!')
        end

        if(tosave) 
            fname = ['data/', var, '/',scenario, '_MM'];
            delete([fname, '.mat']);
            File = matfile(fname, 'Writable', true);
            fprintf("Writing file %s\n", fname);
            File.MMs(1:nMM, 1:s2, 1:s3) = MM.MMs; 
            File.models(1:nMM,1:2) = MM.models; 
            File.trust (1:nMM,1) = MM.trust; 
            NanSims = h.model(any(any(isnan(h.runs),2),3),:);
            File.NanSims = NanSims;
            File.time = h.time(1,:);
            if(isfield(h, 'indices'))
                File.indices = h.indices;
            end
        end

        if(~strcmp(realm, 'amip'))
            piC = load(['data/', var, '/', piCs, '_all.mat']); piC.runs(piC.runs==0)=NaN; piC_lengths = sum(~isnan(piC.runs(:,:,1)), 2);
            [~,s2,s3] = size(piC.runs);
            T_piC = table(piC.model, piC.runs, piC.time, piC_lengths, 'VariableNames', {'model', 'runs', 'time', 'length'});
            
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
            MM.piC_trust = sqrt(histcounts(model_groupings, 1:nMM));
            if(tosave)  
                File.piC_MMs(1:nMM, 1:s2,1:s3) = MM.piC_MMs;
                File.piC_models(1:nMM,1:2) = MM.piC_models; 
                File.piC_trust (1:nMM,1) = MM.piC_trust;
            end
        end

        [~,s2,s3] = size(MM.MMs);
        [GM.models, ~, model_groupings] = unique(MM.models(:,1)); nGM = max(model_groupings);
        weights=splitapply(vert_sum, MM.trust, model_groupings);
        GM.GMs = splitapply(vert_sum, MM.trust.*MM.MMs./weights(model_groupings), model_groupings);
        GM.trust = splitapply(@sum, MM.trust, model_groupings)./sqrt(histcounts(model_groupings, 1:nGM));
        if(tosave)
            fname = ['data/',var, '/', scenario, '_GM'];
            delete([fname, '.mat']);
            File = matfile(fname, 'Writable', true);
            fprintf("Writing file %s\n", fname);
            File.GMs(1:nGM, 1:s2, 1:s3) = GM.GMs; 
            File.models(1:nGM,1) = GM.models;
            File.trust (1:nGM,1) = GM.trust; 
            File.MMM(1, 1:s2, 1:s3) = GM.MMM;
            File.time=h.time(1,:);
            if(isfield(h, 'indices'))
                File.indices = h.indices;
            end
        end
        if(~strcmp(realm, 'amip'))
            [~,s2,s3] = size(MM.piC_MMs);
            [GM.piC_models, I, model_groupings] = unique(MM.piC_models(:,1)); nGM = max(model_groupings);
            weights=splitapply(vert_sum, MM.piC_trust, model_groupings);
            GM.piC_GMs = splitapply(vert_sum, MM.piC_trust.*MM.piC_MMs./weights(model_groupings), model_groupings);
            GM.piC_trust = splitapply(@sum, MM.piC_trust, model_groupings)./sqrt(histcounts(model_groupings, 1:nGM));
            if(tosave)
                File.piC_GMs(1:nGM,1:s2, 1:s3) = GM.piC_GMs;
                File.piC_models(1:nGM,1) = GM.piC_models;
                File.piC_trust(1:nGM,1) = GM.piC_trust;
            end
        end
    end
end

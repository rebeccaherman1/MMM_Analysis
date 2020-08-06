%saves weighted mean anomalies. (Why do I say it's anomalies? It's not!)
tosave = true;

scenarios = {'cmip6_h','cmip6_a', 'cmip6_n', 'cmip6_g'};%'cmip6_h', 'v'};%'r'};%'a6'};%'e'};%'h','a','n','g'};%'amip'};%, 
%vert_nanmean = @(X) nanmean(X,1); 
vert_mean = @(X) mean(X,1); vert_sum = @(X) sum(X,1);
start_year = 1901;%50;

for j = 1:length(scenarios)
    scenario = scenarios{j};
    fprintf("Accessing scenario %s\n", scenario);
    h = load(['data/', scenario, '_all.mat']); s = size(h.runs);
    [model_names, I, model_groupings] = unique(h.model(:,2)); nMM = max(model_groupings);
    MM.MMs = splitapply(vert_mean, h.runs, model_groupings);
    MM.models = [h.model(I,1), model_names];
    MM.trust = sqrt(histc(model_groupings(~any(isnan(h.runs),2)), 1:nMM));
    
    if(tosave) 
        fname = ['data/',scenario, '_MM'];
        File = matfile(fname, 'Writable', true);
        fprintf("Writing file %s\n", fname);
        File.MMs(1:nMM, 1:s(2)) = MM.MMs; 
        File.models(1:nMM,1:2) = MM.models; 
        File.trust (1:nMM,1) = MM.trust; 
        NanSims = h.model(any(isnan(h.runs),2),:)
        File.NanSims = NanSims;
    end
    
    %use for historical simulations; comment out for amip simulations. 
    piC = load(['data/cmip6_piC_all.mat']); piC.runs(piC.runs==0)=NaN; s = size(piC.runs);
    relevant_pC_models = ismember(piC.model(:,1),h.model(:,1));
    piC.model = piC.model(relevant_pC_models,:); piC.runs = piC.runs(relevant_pC_models, :);
    [model_names, I, model_groupings] = unique(piC.model(:,2)); nMM = max(model_groupings);
    MM.piC_MMs = splitapply(vert_mean, piC.runs, model_groupings);
    MM.piC_models = [piC.model(I,1), model_names]; 
    MM.piC_trust = sqrt(histc(model_groupings, 1:nMM));
    if(tosave)  
        File.piC_MMs(1:nMM, 1:s(2)) = MM.piC_MMs;
        File.piC_models(1:nMM,1:2) = MM.piC_models; 
        File.piC_trust (1:nMM,1) = MM.piC_trust;
    end

    s = size(MM.MMs);
    [GM.models, ~, model_groupings] = unique(MM.models(:,1)); nGM = max(model_groupings);
    weights=splitapply(vert_sum, MM.trust, model_groupings);
    GM.GMs = splitapply(vert_sum, MM.trust.*MM.MMs./weights(model_groupings), model_groupings);
    GM.trust = splitapply(@sum, MM.trust, model_groupings)./sqrt(histc(model_groupings, 1:nGM));
    GM.MMM = mean(GM.trust.*GM.GMs/mean(GM.trust), 1);
    if(tosave)
	    fname = ['data/',scenario, '_GM'];
        File = matfile(fname, 'Writable', true);
        fprintf("Writing file %s\n", fname);
        File.GMs(1:nGM, 1:s(2)) = GM.GMs; 
        File.models(1:nGM,1) = GM.models;
        File.trust (1:nGM,1) = GM.trust; 
        File.MMM(1, 1:s(2)) = GM.MMM;
    end	  
    s = size(MM.piC_MMs);
    [GM.piC_models, I, model_groupings] = unique(MM.piC_models(:,1)); nGM = max(model_groupings);
    weights=splitapply(vert_sum, MM.piC_trust, model_groupings);
    GM.piC_GMs = splitapply(vert_sum, MM.piC_trust.*MM.piC_MMs./weights(model_groupings), model_groupings);
    GM.piC_trust = splitapply(@sum, MM.piC_trust, model_groupings)./sqrt(histc(model_groupings, 1:nGM));
    if(tosave)
        File.piC_GMs(1:nGM,1:s(2)) = GM.piC_GMs;
        File.piC_models(1:nGM,1) = GM.piC_models;
        File.piC_trust(1:nGM,1) = GM.piC_trust;
    end
end

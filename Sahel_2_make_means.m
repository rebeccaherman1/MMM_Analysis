%saves weighted mean anomalies. (Why do I say it's anomalies? It's not!)
tosave = true;

scenarios = {'cmip6_h'};%'v'};%'r'};%'a6'};%'e'};%'h','a','n','g'};%'amip'};%, 
vert_mean = @(X) nanmean(X,1); vert_sum = @(X) nansum(X,1);
start_year = 1901;%50;

for j = 1:length(scenarios)
    scenario = scenarios{j};
    fprintf("Accessing scenario %s\n", scenario);
    h = load(['data/', scenario, '_all.mat']); s = size(h.runs);
    [MM.model_names, I, model_groupings] = unique(h.model(:,2)); nMM = max(model_groupings);
    MM.MMs = splitapply(vert_mean, h.runs, model_groupings);
    MM.models = h.model(I,1);
    MM.trust = sqrt(histc(model_groupings(~isnan(h.runs(:,1))), 1:nMM));
    
    if(tosave) 
        fname = ['data/',scenario, '_MM'];
        File = matfile(fname, 'Writable', true);
        fprintf("Writing file %s\n", fname);
        File.MMs(1:nMM, 1:s(2)) = MM.MMs; 
        File.models(1:nMM,2) = MM.model_names; File.models(1:nMM,1) = MM.models; 
        File.trust (1:nMM,1) = MM.trust; 
    end
%{
    piC = load(['data/piC_all.mat']); piC.runs(piC.runs==0)=NaN; s = size(piC.runs);
    relevant_pC_models = ismember(piC.model(:,1),h.model(:,1));
    piC.model = piC.model(relevant_pC_models,:); piC.runs = piC.runs(relevant_pC_models, :);

    [model_names, I, model_groupings] = unique(piC.model(:,2)); nMM = max(model_groupings);
    File.piC_MM(1:nMM, 1:s(2)) = splitapply(vert_mean, piC.runs, model_groupings);
    File.piC_models(1:nMM,2) = model_names; File.piC_models(1:nMM,1) = piC.model(I,1);
    File.piC_trust (1:nMM,1) = sqrt(histc(model_groupings, 1:nMM));
%}

    s = size(MM.MMs);
    [GM.model_names, ~, model_groupings] = unique(MM.models(:,1)); nGM = max(model_groupings);
    weights=splitapply(vert_sum, MM.trust, model_groupings);
    GM.GMs = splitapply(vert_sum, MM.trust.*MM.MMs./weights(model_groupings), model_groupings);
    GM.trust = splitapply(@sum, MM.trust, model_groupings)./sqrt(histc(model_groupings, 1:nGM));
    MMM = mean(GM.trust.*GM.GMs/mean(GM.trust), 1);

    if(tosave)
	fname = ['data/',scenario, '_GM'];
        File = matfile(fname, 'Writable', true);
        fprintf("Writing file %s\n", fname);
        File.GMs(1:nGM, 1:s(2)) = GM.GMs; 
        File.models(1:nGM,1) = GM.model_names;
        File.trust (1:nGM,1) = GM.trust; 
        File.MMM(1, 1:s(2)) = MMM;
    end	
end

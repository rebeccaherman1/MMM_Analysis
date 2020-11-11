cmip5_scenarios = {'h', 'a', 'n', 'g'};
cmip6_scenarios = {'cmip6_h', 'cmip6_a', 'cmip6_n', 'cmip6_g'};
amip_scenarios = {'amip-hist', 'amip-piF'};

r = cmip6_scenarios;
L = length(r);
sims = cell(1, L);
for s = 1:L
    sims(s) = {load(['data/ts/', r{s}, '_all.mat'])};
end
if(L>=4)
    common_models = intersect(intersect(unique(sims{1}.model(:,1)), unique(sims{2}.model(:,1))), unique(sims{4}.model(:,1))); 
else
    common_models = union(unique(sims{1}.model(:,1)), unique(sims{2}.model(:,1)));
end

Tbls = cell(1,L);
models_holder = cell(1,L);
%[Insts_and_Models, ~, master_groupings] = unique(sims{1}.model(ismember(sims{1}.model(:,1), common_models),1:2)); %ALL simulations are inclusive.

for s = 1:L
    to_use = ismember(sims{s}.model(:,1), common_models);
    M = sims{s}.model(to_use,:);
    [mdls, ia, groupings] = unique(M(:,2));
    insts = M(ia,1);
    num_runs = splitapply(@length, groupings, groupings);
    T = table(insts, mdls, num_runs, 'VariableNames', {'Institutions', 'Models', [strrep(r{s}, '-', '_'), '_runs']});
    Tbls(s) = {T};
end

T = Tbls{1};
for i = 2:L
    T = outerjoin(T, Tbls{i}, 'Keys', {'Institutions', 'Models'}, 'MergeKeys', true);
end

T(:,3:end) = fillmissing(T(:,3:end), 'constant', 0)
Num_Runs = sum(T{:,3:end},1)
Num_Models = sum(T{:,3:end}~=0,1)
[~,~,inst_grpngs] = unique(T.Institutions);
Num_Insts = sum(splitapply(@(x) any(x,1), T{:,3:end}~=0, inst_grpngs), 1)


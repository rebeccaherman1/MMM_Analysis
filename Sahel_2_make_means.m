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

global h_indices, global lat, global lon

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

for j = 1:length(scenarios)
    if(gts)
    	AA = load(make_data_filename('pr', 7,9, CM, 'all'));%'MM'));%
    else
        AA = mk_tbl_joined(variables, start_month, end_month, CM, 'all');
    end
    %AA = load(make_data_filename('ts', start_month, end_month, 'cmip6_h', 'all'));%'a', 'MM'));%
    common_models = unique(AA.institution); %just model for CMIP6?
    
    clear MM GM
    scenario = scenarios{j};
    fprintf("Accessing scenario %s\n", scenario);% variable %s, var);

    h = mk_tbl_joined(variables, start_month, end_month, scenario, 'all');
    h = h(ismember(h.institution, common_models),:);
    [MM, model_names, num_runs] = mk_mns(h, variables, vert_nan_mean, '_MMs');
    if(any(h.time~=h.time(1,:)))
        fprintf('time does not line up!')
    end

    if(~strcmp(realm, 'amip') && ~contains(var, 'global'))
        piC = mk_tbl_joined(variables, start_month, end_month, piCs, 'all');
        if ismember('time', piC.Properties.VariableNames)
            piC = removevars(piC, 'time');
        end
        piC_lengths = nan(height(piC), length(variables));
        for v = 1:length(variables)
            var = variables{v};
            piC.(var)(piC.(var)==0)=NaN; piC_lengths(:,v) = sum(~isnan(piC.(var)(:,:,1)), 2);
        end
        lp = sum(any(piC_lengths~=piC_lengths(:,1),2));
        if lp>0
            fprintf("length of piC runs differs between variables for %i simulations\n", lp)
        end
        piC_lengths = min(piC_lengths, [], 2);
        piC.length = piC_lengths;
        
        %hacky name inconsistency changes; should actually remove after
        %re-downloading everything. But I'm still not able to re-download
        %the piC simulations. Which... means they could be wrong... 
        if(any(strcmp(piC.institution, 'CESM')))
            piC.institution(strcmp(piC.institution, 'CESM'),1)={'NCAR'};
        end
        if(any(strcmp(piC.institution, 'CanESM')))
            piC.institution(strcmp(piC.institution, 'CanESM'),1)={'CCCma'};
        end
        if(any(strcmp(piC.institution, 'GISS')))
            piC.institution(strcmp(piC.institution, 'GISS'),1)={'NASA'};
        end
        if(any(strcmp(piC.institution, 'NorESM')))
            piC.institution(strcmp(piC.institution, 'NorESM'),1)={'Nor'};
        end

        %Use common institutions for CMIP5 but common models for CMIP6
        if(strcmp(realm, 'cmip5'))
            npcm = 'institution'; M_N = unique(MM.model(:,1));
        else
            npcm = 'model'; M_N = model_names;
        end
        relevant_pC_models = ismember(piC.(npcm),M_N);
        piC = piC(relevant_pC_models, :); 
        %before I do this! let me pretend I have more runs!!!!!
        [Lia, Locb] = ismember(piC.model, model_names);
        Locb2 = Locb; Locb2(~Lia)=1; %momentarily lose the 0s so that matlab doesn't crash
        num_pC_reps = num_runs(Locb2); num_pC_reps(~Lia)=1; %replace those numbers with just one repetition
        piC = repelem(piC, num_pC_reps, 1);
        offsets = arrayfun(@randi, piC.length);
        %now I will scramble them a bit!
        for row = 1:size(piC.model, 1)
            lr = piC.length(row); or = offsets(row);
            if(or>1)
                for v = 1:length(variables)
                    var = variables{v};
                    piC.(var)(row, 1:lr) = [piC.(var)(row, or:lr), piC.(var)(row, 1:or-1)];
                end
            end
            hts = size(h.time,2);
            if(lr<hts)
                fprintf('Extending piC simulation %s %s\n',piC.model{row}, piC.run{row})
                for v = 1:length(variables)
                    var = variables{v};
                    piC.(var)(row, (lr+1):hts) = piC.(var)(row, 1:(hts-lr));
                end
            end
        end
        MM = mk_mns(piC, variables, vert_mean, '_piC_MMs', MM);
        flds = fieldnames(MM); piC_flds = contains(flds, 'piC');
        T_piC = mk_tbl(rmfield(MM, flds(~piC_flds)));

        %remove clim values. No need to take anomalies of sims ever
        %again!
        [has_piC, usd] = ismember(MM.model(:,2), MM.piC_model(:,2));
        for v = 1:length(variables)
            var = variables{v};
            MM.([var,'_clim'])(has_piC,:,:) = nanmean(MM.([var,'_piC_MMs'])(usd(has_piC), :), 2);
        end
        %for models with no piC simulation, estimate clim using an
        %average.
        if(any(~has_piC))
            leftovers = MM.model(~has_piC,1);
            lo = nan(length(leftovers), length(variables));
            for v = 1:length(variables)
                var = variables{v};
                for l = 1:length(leftovers)
                    lo(l,v) = nanmean(MM.([var,'_piC_MMs'])(strcmp(MM.piC_model(:,1), leftovers{l}),:),[1,2]);
                end
                MM.([var,'_clim'])(~has_piC,:,:) = lo(:,v);%cell2mat(leftovers);
                MM.([var,'_MMs']) = MM.([var,'_MMs']) - MM.([var,'_clim']);
            end
        end
    end
    
    if(tosave) 
        %TODO update file naming because no longer by variable
        fname = make_data_filename_all(start_month, end_month, scenario, 'MM');
        delete(fname);
        fprintf("Writing file %s\n", fname);
        %MM.NanSims = h.model(any(any(isnan(h.runs),2),3),:);
        if ismember('time', h.Properties.VariableNames)
            MM.time = h.time(1,:);
        end
        if(exist('h_indices', 'var'))
            MM.indices = h_indices;
        end
        if(isfield('h', 'lat'))
            MM.lat = h.lat;
            MM.lon = h.lon;
        end
        save(fname, '-struct', 'MM')
    end

    MM_T = mk_tbl(rmfield(MM, flds(piC_flds)));
    GM = mk_g_mns(MM_T, variables,vert_sum,'_GMs');
    if(~strcmp(realm, 'amip') && ~any(contains(variables, 'global')))
        GM = mk_g_mns(T_piC, variables,vert_sum,'_piC_GMs', GM);
    end
    if(tosave)
        fname = make_data_filename_all(start_month, end_month, scenario, 'GM');
        delete(fname);
        fprintf("Writing file %s\n", fname);
        GM.time=h.time(1,:);
        if(exist('h_indices', 'var'))
            GM.indices = h_indices;
        end
        if(sum(size(lat))>0)
            GM.lat = lat;
            GM.lon = lon;
        end
        save(fname, '-struct', 'GM')
    end
end

function [ht] = mk_tbl(h)
    global h_indices, global lat, global lon
    flds = fieldnames(h);
    M = flds{contains(flds, 'model')};
    ht = table(h.(M)(:,1), h.(M)(:,2),... h.model(:,3),h.runs,...
        'VariableNames', {'institution', 'model'});%, 'run', 'runs'});
    if(size(h.(M),2)>2)
        ht.run = h.(M)(:,3);
    end
    for f = 1:length(flds)
        fld = flds{f};
        if(strcmp(fld, M))
        elseif(strcmp(fld, 'time'))
            ht.time = repmat(h.time, size(h.model, 1)/size(h.time,1),1);
        elseif(strcmp(fld, 'indices'))
            h_indices = h.indices(1,:);
        elseif(strcmp(fld, 'lat'))
            lat = h.lat;
            lon = h.lon;
        else
            ht.(fld) = h.(fld);
        end
    end
end

%TODO combine to single source of truth
%unique input, nm of target var
function [GM, model_names] = mk_g_mns(MM, variables,f,nm,GM)
    [model_names, I, model_groupings] = unique(MM.institution); nMM = max(model_groupings);
    if(nargin < 5)
        nm1 = 'model';
        nm2 = 'trust';
        do_MMM = true;
    else
        nm1 = 'piC_model';
        nm2 = 'piC_trust';
        do_MMM = false;
    end
    weights=splitapply(f, MM.(nm2), model_groupings);
    GM.(nm1) = MM.institution(I);
    GM.(nm2) = splitapply(@sum, MM.(nm2), model_groupings)./sqrt(histcounts(model_groupings, (0:nMM)+.5)');
    for v = 1:length(variables)
        vbls = MM.Properties.VariableNames;
        var = variables{v}; 
        v_MM = vbls{contains(vbls, var) & contains(vbls, 'MM')};
        GM.([var,nm]) = splitapply(f, MM.(nm2).*MM.(v_MM)./weights(model_groupings), model_groupings);
        if(do_MMM)
            GM.([var,'_clim']) = splitapply(f, MM.(nm2).*MM.([var,'_clim'])./weights(model_groupings), model_groupings);
            GM.([var,'_MMM']) = sum(GM.(nm2).*GM.([var,nm]), 1)/sum(GM.(nm2));
        end
    end
end

function [MM, model_names, num_runs] = mk_mns(h, variables,f,nm,MM)
    [model_names, I, model_groupings] = unique(h.model); nMM = max(model_groupings);
    some_non_nans = varfun(@(X) ~all(all(isnan(X),2),3), h(:,ismember(h.Properties.VariableNames, variables)));
    num_runs = histcounts(model_groupings(all(some_non_nans{:,:},2)), (0:nMM)+.5)';
    for v = 1:length(variables)
        var = variables{v};
        MM.([var,nm]) = splitapply(f, h.(var), model_groupings);
    end
    if(nargin < 5)
        nm1 = 'model';
        nm2 = 'trust';
    else
        nm1 = 'piC_model';
        nm2 = 'piC_trust';
    end
    MM.(nm1) = [h.institution(I), model_names];
    MM.(nm2) = sqrt(num_runs);
end

function [AA] = mk_tbl_joined(vars, start_month, end_month, CM, type)
    L = length(vars); K = {'institution', 'model', 'run'};
    D = cell(L,1); 
    for i = 1:L
        D{i} = mk_tbl(load(make_data_filename(vars{i}, start_month, end_month, CM, type)));%'a', 'MM'));%
    end
    AA = D{1}; AA = update_var_name(AA,'runs',vars{1});
    te = ismember('time', AA.Properties.VariableNames);
    if te
        AA = update_var_name(AA, 'time', 'time_AA');
    end
    for i=2:L
        AA = innerjoin(AA, D{i}, 'Keys', K);
        if ismember('time', AA.Properties.VariableNames)
            if ~all(AA.time_AA==AA.time, 'all')
                fprintf("Time Mismatch!\n")
            end
            AA = removevars(AA, {'time'});
        else
            fprintf("Some variables don't have time\n")
        end
        AA = update_var_name(AA, 'runs', vars{i});
    end
    if te
        AA = update_var_name(AA, 'time_AA', 'time');
    end
end

function [T] = update_var_name(T, v_old, v_new)
    T.Properties.VariableNames{strcmp(T.Properties.VariableNames, v_old)}=v_new;
end

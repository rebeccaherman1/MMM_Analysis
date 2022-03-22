%TODO fix this script so that it saves ALL VARIABLES that have the right size...!

scenarios = {'cmip6_piC'};%, 'cmip6_a', 'cmip6_n', 'cmip6_g'};%cle'amip-hist'};%
remove_models = {...
    %'E3SM-1-1 p1'};%, ... Inf, model name
    %'MIROC-ES2H', 'KIOST-ESM p1'};% zg infinity or 0. model name.
    'MIROC-ES2H'};% ta  0. model name.
    %'MIROC-ES2', 'GISS-E2'};%hus_bndries, infty and 0
    %'MIROC-ES2H', 'GISS-E2'};%hus, infty and 0
    %'NASA-GISS'};%,...
    %'MOHC'};% 0, inst name
    % {... prcp values 3 orders of magnitude too...
    %'CIESM',... small, in piC and historical simulations. 
    %'MCM-UA-1-0'}; %large, in piC simlations.
variables = {'ta'};%'M_va', 'penetration'};%_conv_925'};%'pr', 
start_month = 7;
end_month = 9;
for v = 1:length(variables)
    for s = 1:length(scenarios)
        fname = make_data_filename(variables{v}, start_month, end_month, scenarios{s}, 'all');
        h = load(fname); 
        to_remove = zeros(size(h.model,1),1);
        if(istable(h.model))
            M = h.model{:,1:2};
        else
            M = h.model(:,1:2);
        end
        for rm = 1:length(remove_models)
            to_remove = to_remove | any(contains(M, remove_models{rm}), 2);
        end
        fprintf([scenarios{s}, ', ', variables{v}, ': removing models:'])
        h.model(to_remove,:)
        for fld = fieldnames(h)'
            if(istable(h.(fld{:})))
                h.(fld{:}) = h.(fld{:})(~to_remove,:);
            else
                h.(fld{:}) = h.(fld{:})(~to_remove,:,:);
            end
        end
        save(fname, '-struct', 'h');
        %{
        model = h.model(~to_remove,:);
        runs = h.runs(~to_remove,:,:);
        time = h.time(~to_remove,:);
        if(isfield(h, 'indices'))
            indices = h.indices;
            save(fname, 'model', 'runs', 'time', 'indices')
        else
            save(fname,'model','runs','time');
        end
        %}
    end
end

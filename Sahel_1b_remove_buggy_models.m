scenarios = {'cmip6_h'};%, 'cmip6_a', 'cmip6_n', 'cmip6_g'};%cle'amip-hist'};%
remove_models = ...
    {... hus_conv_925
    'E3SM-1-1 p1' ... Inf, model name
    'MOHC'};% 0, inst name
    % {... prcp values 3 orders of magnitude too...
    %'CIESM',... small, in piC and historical simulations. 
    %'MCM-UA-1-0'}; %large, in piC simlations.
variables = {'hus_conv_925'};%'pr', 
for v = 1:length(variables)
    for s = 1:length(scenarios)
        fname = make_model_filename(variables{v}, start_month, end_month, scenarios{s}, 'all');
        h = load(fname); sz = size(h.runs);
        to_remove = zeros(sz(1),1);
        for rm = 1:length(remove_models)
            to_remove = to_remove | any(contains(h.model(:,1:2), remove_models{rm}), 2);
        end
        fprintf([scenarios{s}, ': removing models:'])
        h.model(to_remove,:)
        model = h.model(~to_remove,:);
        runs = h.runs(~to_remove,:,:);
        time = h.time(~to_remove,:);
        if(isfield(h, 'indices'))
            indices = h.indices;
            save(fname, 'model', 'runs', 'time', 'indices')
        else
            save(fname,'model','runs','time');
        end
    end
end

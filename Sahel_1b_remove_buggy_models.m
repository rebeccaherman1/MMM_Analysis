scenarios = {'cmip6_h', 'cmip6_a', 'cmip6_n', 'cmip6_g'};%cle'amip-hist'};%
remove_models = {... prcp values 3 orders of magnitude too...
    'CIESM',... small, in piC and historical simulations. 
    'MCM-UA-1-0'}; %large, in piC simlations.
variables = {'pr', 'ts'};
for v = 1:length(variables)
    for s = 1:length(scenarios)
        fname = ['data/', variables{v}, '/', scenarios{s}, '_all.mat'];
        h = load(fname); sz = size(h.runs);
        to_remove = zeros(sz(1),1);
        for rm = 1:length(remove_models)
            to_remove = to_remove | contains(h.model(:,2), remove_models{rm});
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

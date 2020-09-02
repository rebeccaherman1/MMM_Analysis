scenarios = {'cmip6_r'};%'cmip6_h', 'cmip6_a', 'cmip6_n', 'cmip6_g'};
remove_models = {'CIESM', 'MCM-UA-1-0'};
for s = 1:length(scenarios)
    fname = ['data/pr/', scenarios{s}, '_all.mat'];
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
    save(fname,'model','runs','time');
end

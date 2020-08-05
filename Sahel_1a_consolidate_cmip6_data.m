clear
model_file_name = 'data/cmip6_h_all.mat';
folder = '~/netcdf/cmip6/historical';
files = split(ls(folder));
next_line = 1;
for file = files(1:end-1)'
    names = make_model_and_p_name(file);
    model(next_line, 1:3) = names;
    fopen_name = [folder, '/', file{:}];
    pr = ncread(fopen_name, 'pr');
    l = length(pr);
    runs(next_line,1:l)   = pr;
    time(next_line,1:l)   = ncread(fopen_name, 'year');
    next_line=next_line+1;
    save(model_file_name,'model','runs', 'time');
end

function [model_names] = make_model_and_p_name(file)
% special version for PSL-FACTS
    fields = split(file, [".", "_"]);
    institution = fields{2};
    model = fields{3};
    run = fields{4};
    run_stats = split(run, ["p", "f"]);
    model_and_p = [model, ' p', run_stats{2}];
    model_names = {institution, model_and_p, run};
end

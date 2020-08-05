clear
model_file_name = 'data/cmip6_h_all.mat';
folder = '~/netcdf/cmip6/historical';
files = split(ls(folder));
next_line = 1;
T = 1850:1:2014;
for file = files(1:end-1)'
    names = make_model_and_p_name(file);
    model(next_line, 1:3) = names;
    fopen_name = [folder, '/', file{:}];
    Time = ncread(fopen_name, 'year');
    l = length(T);
    if Time(1)<=1901 && Time(end)>=2003
	T_x = ismember(Time, T);
	pr = ncread(fopen_name, 'pr');
        runs(next_line,1:l)   = pr(T_x); 
        time(next_line,1:l)   = T(T_x);
        next_line=next_line+1;
        save(model_file_name,'model','runs', 'time');
    end
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

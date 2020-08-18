clear
variable = 'ts';

scenarios = {'historical'};%, 'hist-aer', 'hist-nat', 'hist-GHG', 'piControl'};
short_names = {'h', 'a', 'n', 'g', 'piC'};
for i = 1:length(scenarios)
    model_file_name = ['data/', variable, '/cmip6_', short_names{i}, '_all.mat']
    folder = ['~/netcdf/cmip6/preprocessed/', scenarios{i}];
    files = split(ls(folder));
    files = files(contains(files, [variable, '_']));
    next_line = 1;
    T = 1850:1:2014;
    for file = files'
        names = make_model_and_p_name(file);
        model(next_line, 1:3) = names;
        fopen_name = [folder, '/', file{:}];
        %ncdisp(fopen_name)
        INFO = ncinfo(fopen_name);
        if any(contains({INFO.Variables.Name}, {'time'}))
        	Time = ncread(fopen_name, 'time');
        else
            Time = ncread(fopen_name, 'year');
        end
        %pr = ncread(fopen_name, 'pr');
        NA = ncread(fopen_name, 'NA');
        GT = ncread(fopen_name, 'GT');
        NARI = ncread(fopen_name, 'NARI');
        pr = cat(3, NA, GT, NARI);
        if contains(model_file_name, 'piC')
            l = length(Time);
        elseif Time(1)<=1901 && Time(end)>=2014
            T_x = (Time >= 1901) & (Time <= 2014);
            Time = Time(T_x);
            pr = pr(T_x,:,:);
            l = length(pr);
            if(length(Time)~=length(pr))
                fprintf('pr contains fill values');
            end
        else
            continue
        end
        runs(next_line,1:l,1:3)   = pr; 
        time(next_line,1:l)   = Time;
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

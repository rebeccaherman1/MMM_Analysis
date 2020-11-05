%Only use on PREPROCESSED data which was downloaded from the cloud via my
%jupyter script. For data I downloaded in full, see matlab code
%Sahel_1_save_data_amip.

clear
variable = 'pr';

scenarios = {'historical'};%'amip-hist', 'hist-aer', 'hist-nat', 'hist-GHG', 'piControl'};
short_names = {'cmip6_h'};%'amip-hist', 'cmip6_a', 'cmip6_n', 'cmip6_g', 'cmip6_piC'};

for i = 1:length(scenarios)
    clear model runs time
    model_file_name = ['data/', variable, '/', short_names{i}, '_all.mat']
    folder = ['~/netcdf/cmip6/preprocessed/', scenarios{i}];
    files = split(ls(folder));
    files = files(contains(files, [variable, '_']));
    next_line = 1;
    for file = files'
        names = make_model_and_p_name(file);
        model(next_line, 1:3) = names;
        fopen_name = [folder, '/', file{:}];
        %ncdisp(fopen_name)
        INFO = ncinfo(fopen_name);
        if any(contains({INFO.Variables.Name}, {'time'}))
            Time = ncread(fopen_name, 'time');
        elseif any(contains({INFO.Variables.Name}, {'year'}))
            Time = ncread(fopen_name, 'year');
        end
        if(strcmp(variable, 'pr'))
            pr = ncread(fopen_name, 'pr');
            s3=1;
        else
            NA = ncread(fopen_name, 'NA');
            GT = ncread(fopen_name, 'GT');
            NARI = ncread(fopen_name, 'NARI');
            pr = cat(3, NA, GT, NARI);
            s3=3;
            indices = permute({'NA', 'GT', 'NARI'},[1,3,2]);
        end
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
        runs(next_line,1:l,1:s3) = pr; 
        time(next_line,1:l) = Time;
        next_line=next_line+1;
        if(strcmp(variable, 'ts'))
            save(model_file_name, 'model', 'runs', 'time', 'indices')
        else
            save(model_file_name,'model','runs', 'time');
        end
    end
end

function [model_names] = make_model_and_p_name(file)
% special version for PSL-FACTS
    fields = split(file, [".", "_"]);
    model = fields{3};
    run = fields{4};
    run_stats = split(run, ["p", "f"]);
    model_and_p = [model, ' p', run_stats{2}];
    institution = fields{2};
    model_names = {institution, model_and_p, run};
end

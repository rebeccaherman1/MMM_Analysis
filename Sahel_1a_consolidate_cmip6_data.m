%Only use on PREPROCESSED data which was downloaded from the cloud via my
%jupyter script. For data I downloaded in full, see matlab code
%Sahel_1_save_data_amip.

%TODO check why CSIRO files give segfaults.

clear
variable = 'ta';
location = 'Sahel';

scenarios = {'historical'};%, 'hist-aer'};%, 'hist-nat', 'hist-GHG'};%'amip-hist', 'piControl'};
short_names = {'cmip6_h', 'cmip6_a', 'cmip6_n', 'cmip6_g'};%'amip-hist', 'cmip6_piC'};

for i = 1:length(scenarios)
    clear model runs time
    fldr_name = ['data/', variable];
    if strcmp(variable, 'zg') || strcmp(variable, 'ta')
	fldr_name = [fldr_name, '100'];
    end
    if ~exist(fldr_name, 'dir')
	mkdir(fldr_name);
    end
    model_file_name = [fldr_name, '/', short_names{i}, '_all.mat']
    folder = ['~/netcdf/cmip6/preprocessed/', scenarios{i}];
    files = split(ls(folder));
    files = files(contains(files, [variable, '_']));
    next_line = 1;
    for file = files'
        names = make_model_and_p_name(file, variable);
        model(next_line, 1:3) = names;
        fopen_name = [folder, '/', file{:}];
        %ncdisp(fopen_name)
        INFO = ncinfo(fopen_name);
	if strcmp(variable, 'ts') && ~any(contains({INFO.Variables.Name}, {'SA'})) || contains(fopen_name, 'CSIRO')
	    fprintf("skipping file %s\n", file{:})
	    continue;
        else
	    fprintf("processing file %s\n", file{:})
	end
        if any(contains({INFO.Variables.Name}, {'time'}))
            Time = ncread(fopen_name, 'time');
        elseif any(contains({INFO.Variables.Name}, {'year'}))
            Time = ncread(fopen_name, 'year');
        end
        if(strcmp(location, 'Sahel'))
	    fprintf('fopen_name = %s, variable = %s\n', fopen_name, variable)
	    V = split(variable, ["_"]);
	    if(any(contains(V, 'conv')))
		V = '__xarray_dataarray_variable__';
	    else
		V = V{1};
	    end
	    pr = ncread(fopen_name, V);
	    s3=1;
        else
            NA = ncread(fopen_name, 'NA');
            GT = ncread(fopen_name, 'GT');
            NARI = ncread(fopen_name, 'NARI');
            p1 = NA + GT;
            SA = ncread(fopen_name, 'SA');
            md = ncread(fopen_name, 'md');
            pr = cat(3, NA', GT', NARI', SA', md');
            s3=size(pr, 3);
            indices = permute({'NA', 'GT', 'NARI', 'SA', 'md'},[1,3,2]);
        end
        if contains(model_file_name, 'piC')
            l = length(Time);
        elseif Time(1)<=1901 && Time(end)>=2014
            T_x = (Time >= 1901) & (Time <= 2014);
            Time = Time(T_x);
            if(ndims(pr)>2)
		pr = pr(:,T_x,:);
	    elseif (size(pr,2)==1)
	        pr = pr(T_x)';
	    else
		pr = pr(end,T_x); %picking plev = 100 hPa for those files where I took too much. 
	    end
            l = size(pr,2);
            if(length(Time)~=l)
                fprintf('pr contains fill values\n');
            end
        else
            continue
        end
        runs(next_line,1:l,1:s3) = pr; 
        time(next_line,1:l) = Time;
        next_line=next_line+1;
        if(strcmp(location, 'Ocean'))
            save(model_file_name, 'model', 'runs', 'time', 'indices')
        else
            save(model_file_name,'model','runs', 'time');
        end
    end
end

%{
CMIP6 institutions file calculated using the following code, after 
downloading CMIP6 historical simulations using jupyter:
%}
%{
H = load('data/pr/cmip6_h_all.mat');
[~, ia, ~] = unique(H.model(:,2));
institutions = H.model(ia, 1:2);
save('data/institutions_cmip6.mat', 'institutions');
%}

function [model_names] = make_model_and_p_name(file, var)
% special version for PSL-FACTS
    var_length = length(split(var, [".", "_"]));
    fields = split(file, [".", "_"]);
    model = fields{2+var_length};
    run = fields{3+var_length};
    run_stats = split(run, ["p", "f"]);
    model_and_p = [model, ' p', run_stats{2}];
    institution = fields{1+var_length};
    model_names = {institution, model_and_p, run};
end

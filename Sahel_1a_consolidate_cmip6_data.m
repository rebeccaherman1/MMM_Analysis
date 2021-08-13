%Only use on PREPROCESSED data which was downloaded from the cloud via my
%jupyter script. For data I downloaded in full, see matlab code
%Sahel_1_save_data_amip.
%variables that are included in some models but not others can be found in SKIPPED_VARS
%TODO check why CSIRO files give segfaults.

clear
variable = 'ts';
%location = 'Sahel'; Not currently used. Perhaps use ~strcmp(location, Sahel) for the ocean basins instead of strcmp(variable, ts)
start_month = 7;
end_month = 9;

scenarios = {'historical','hist-aer','hist-nat', 'hist-GHG','piControl'};%};%};%, , 'amip-hist', 
short_names = {'cmip6_h','cmip6_a','cmip6_n', 'cmip6_g','cmip6_piC'};%};%};%, , 'amip-hist', 
skipped_vars = cell(1,6);

for i = 1:length(scenarios)
    fprintf('\nSCENARIO %s\n', scenarios{i}')
    %create folder and filename where compiled data will be saved
    [model_file_name, fldr_name] = make_data_filename(variable, start_month, end_month, short_names{i}, 'all');
    if ~exist(fldr_name, 'dir')
	mkdir(fldr_name);
    end
    %get list of files to compile
    folder = ['~/netcdf/cmip6/preprocessed/', scenarios{i}, '/', num2str(start_month), '-', num2str(end_month)];
    files = split(ls(folder));
    files = files(contains(files, [variable, '_']));
    for file = files'
        fopen_name = [folder, '/', file{:}];
        %ncdisp(fopen_name)
        INFO = ncinfo(fopen_name);
	vars = {INFO.Variables.Name};
%	if strcmp(variable, 'ts') && ~any(contains(vars, {'SA'})) || contains(fopen_name, 'CSIRO')
%	    fprintf("skipping file %s\n", file{:})
%	    continue;
%        else
	    fprintf("processing file %s\n", file{:})
%	end
	%downloads all available data into a cell array of varying length
	Dims = {INFO.Variables.Dimensions};
	D = cell(length(vars), 1);
	for v = 1:length(vars)
	    d = ncread(fopen_name, vars{v});
	    Lengths = [Dims{v}.Length];
	    dims = {Dims{v}.Name};
	    %make sure the first dimension has length 1 and that time is the second dimension.
	    D{v} = permute(d, [length(Lengths)+1, find(strcmp(dims, 'time') | strcmp(dims, 'year')),...
		                   find(~(strcmp(dims, 'time') | strcmp(dims, 'year')))]);
        end
	%rename 'year' to 'time' if needed
	if(any(contains(vars, 'year')))
	    vars{strcmp(vars, 'year')} = 'time';
        end
	Time = D{strcmp(vars, 'time')};
	%if the variable is ts, combine the basins into one array called VARIABLE.
        if(strcmp(variable, 'ts'))
	    %Add p1 if missing.
	    if(~any(contains(vars, 'p1')) & contains(vars, 'NA') & contains(vars, 'GT'))
		p1 = D{strcmp(vars, 'NA')} + D{strcmp(vars, 'GT')};
                D = [D; {p1}];
	        vars = [vars, {'p1'}];
	    end
	    %keep this list of all potential ocean basins up to date!
	    all_basins = contains(vars, {'NA', 'GT', 'NARI', 'p1', 'md', 'SA', 'TA', 'GG'});
	    pr = cat(3, D{all_basins});
	    indices = vars(all_basins);
	    D = {pr; Time; indices};
	    vars = {variable, 'time', 'indices'};
        end
	%Find the variable of interest
	if(contains(variable, 'conv'))
	    V = '__xarray_dataarray_variable__';
        elseif(contains(variable, 'bndries'))
	    V = variable(1:(end-8));
        else	    
	    V = variable;
        end
        %restrict time to the years of interest.
	if contains(model_file_name, 'piC')
            l = length(Time);
        elseif Time(1)<=1901 && Time(end)>=2014
            T_x = (Time >= 1901) & (Time <= 2014);
            D{strcmp(vars, 'time')} = Time(T_x);
	    Dv = cat(1, D{contains(vars, V)});
	    D(contains(vars, V)) = num2cell(Dv(:,T_x,:), [2,3]);
	    l = size(D{find(contains(vars, V), 1)},2);
            if(length(D{strcmp(vars, 'time')})~=l)
                fprintf('pr contains fill values\n');
            end
        else
            fprintf('does not cover desired time period\n');
	    continue
        end
	%add institution, model, run to the cell array D
	D{length(vars)+1} = make_model_and_p_name(file, variable);
	vars{length(vars)+1} = 'model';
	%concatenate with previous files
	if(~exist('D_tot', 'var'))
	    D_tot = D;
	    vars_tot = vars;
        else
	    for di = 1:length(vars)
		dti = strcmp(vars_tot, vars{di});
		if(any(dti))
		    if(contains(model_file_name, 'piC') & ~any(strcmp(vars{dti}, {'indices', 'model'})))
			tmp = D_tot{dti};
			h1 = size(tmp, 1);
			h2 = size(D{di}, 1);
			w1 = size(tmp, 2);
			w2 = size(D{di}, 2);
			d = size(tmp, 3);
			D_tot{dti} = nan(h1+h2, max(w1, w2));
			D_tot{dti}(1:h1, 1:w1,1:d) = tmp;
			D_tot{dti}((h1+1):(h1+h2), 1:w2, 1:d) = D{di};
		    else
		        D_tot{dti} = cat(1, D_tot{dti}, D{di});
		    end
		else
		    fprintf('missing variable %s\n', vars{di});
		    M = D{strcmp(vars, 'model')};
		    skipped_vars = cat(1, skipped_vars, {scenarios{i}, M{:}, vars{di}, D{di}});
	        end
	    end
        end
    end
    %change the name of the variable of interest to 'runs' to make later code easier. 
    %the name of the variable will be preserved in the MODEL_FILE_NAME
    if contains(variable, 'bndries')
        vars_tot(contains(vars_tot, V)) = arrayfun(@(X) ['runs', X{:}((length(V)+1):end)],...
                                                   vars_tot(contains(vars_tot, V)), 'UniformOutput', 0);
    else
    	vars_tot{strcmp(vars_tot, V)} = 'runs';
    end
    S = cell2struct(D_tot, vars_tot);
    save(model_file_name, '-struct', 'S');
    clear D_tot, vars_tot;
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

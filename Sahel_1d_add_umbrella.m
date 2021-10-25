
variable = 'ts';
sm = 7; em = 9;
scenarios = {'g'};
umbrella=load('data/institutions_cmip5.mat');
for s =1:length(scenarios)
    model_file_name = make_data_filename(variable, sm, em, scenarios{s},'all');
    D = load(model_file_name);
    D.model(:,1) = find_umbrella_names(D.model(:,1), umbrella);
    save(model_file_name, '-struct', 'D');
end

function[umbrella_names] = find_umbrella_names(used_models, umbrella)
   umbrella_names = cell(length(used_models),1);
   %case for cmip6
   if(~isfield(umbrella, 'abbrev'))
       [~,Loc] = ismember(used_models,umbrella.institutions(:,2));
       umbrella_names = umbrella.institutions(Loc,1);
   %case for cmip5
   else
       for k = 1:length(umbrella.models)
           umbrella_names(contains(used_models, umbrella.models(k)), 1) = umbrella.abbrev(k);
       end
   end
end
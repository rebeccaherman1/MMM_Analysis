function [fname] = make_analysis_filename(variable, sc,start_year, end_year, N)
    fname = ['Analysis/', variable, '/', sc, '_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.mat'];
end

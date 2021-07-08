function [fname] = make_analysis_filename(variable, start_month, end_month, sc,start_year, end_year, N)
    fname = ['Analysis/', variable, '/', num2str(start_month), '-', num2str(end_month), '/', sc, '_', num2str(start_year), '-', num2str(end_year), '_N', num2str(N), '.mat']
end

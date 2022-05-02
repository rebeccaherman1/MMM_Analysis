function [fname, fldr] = make_data_filename(variable, start_month, end_month, expt, stage)
    if ischar(start_month)
        fldr = ['data/', variable, '/', start_month, '/'];
    else
        fldr = ['data/', variable, '/', num2str(start_month), '-', num2str(end_month), '/'];
    end
    if(strcmp(expt, 'observations'))
	fl = expt;
    else
	fl = [expt, '_', stage];
    end
    fname = [fldr, fl, '.mat'];
end

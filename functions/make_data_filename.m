function [fname, fldr] = make_data_filename(variable, start_month, end_month, expt, stage)
    fldr = ['data/', variable, '/', num2str(start_month), '-', num2str(end_month), '/'];
    if(strcmp(expt, 'observations'))
	fl = expt;
    else
	fl = [expt, '_', stage];
    end
    fname = [fldr, fl, '.mat'];
end

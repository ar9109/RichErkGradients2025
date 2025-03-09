function [] = plot_fit_gem(analysis_mat,gemFolder,savename_suffix)
arguments
    analysis_mat
    gemFolder
    savename_suffix = [];
end
%PLOT_FIT_GEM Summary of this function goes here
%   Detailed explanation goes here
for i =1:numel(analysis_mat)
    
    disp('------------');
    display(analysis_mat(i).name);

    s = analysis_mat(i);

    ampPointrot = s.ampPointrot;
    endPointrot = s.endPointrot;
    fit_obj = s.(['fit_gem_heavi',savename_suffix]);
    gof = s.(['fit_gem_heavi',savename_suffix,'_gof']);

    if ~isempty(fit_obj) && numel(s.binvalue)>3
        fractionGEM = s.fractionGEM;
        binvalue = s.binvalue;
        
        % ktr
        fit_x = linspace(0,endPointrot(1)-ampPointrot(1));
        binned.x = binvalue;
        binned.y = fractionGEM;
    
        f = figure('visible','off');hold on;

        plot_fit(fit_obj,fit_x,gof,binned = binned,visible='off');
        
        ylabel('%GEM+');
        xlabel('x/\mum')
        xlim([-100 3000]);
        ylim([0 1]);
        title(s.name,'Interpreter', 'none');
        set(gca,'XDir','reverse')

        config_plot(f);

        saveas(f,[gemFolder,filesep,s.name,'_gem',savename_suffix,'.png']);
        close(f);
    end
    

end
end


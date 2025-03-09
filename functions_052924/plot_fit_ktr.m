function [] = plot_fit_ktr(analysis_mat,ktrFolder,savename_suffix)
arguments
    analysis_mat
    ktrFolder
    savename_suffix = [];
end
%PLOT_FIT_KTR Summary of this function goes here
%   Detailed explanation goes here




for i =1:numel(analysis_mat)
    
    disp('------------');
    display(analysis_mat(i).name);

    s = analysis_mat(i);
%     ktr = s.ktr;
%     ccrot = s.ccrot;
    ampPointrot = s.ampPointrot;
    endPointrot = s.endPointrot;
    fit_obj = s.(['fit_ktr_linear',savename_suffix]);
%     gof = s.(['fit_ktr_linear',savename_suffix,'_gof']);

    if ~isempty(fit_obj)
        averageKTR = s.averageKTR;
        semKTR = s.semKTR;
        binvalue = s.binvalue;
        
        % ktr
        fit_x = linspace(0,endPointrot(1)-ampPointrot(1));
        binned.x = binvalue;
        binned.y = averageKTR;
        binned.sem = semKTR;
        data.x = s.ccrot(:,1);
        data.y = s.ktr;

        f = figure('visible','off');hold on;

        plot_fit(fit_obj,fit_x,'a',binned=binned,data = data,visible='off');
        
        ylabel('Erk');
        xlabel('x/\mum')
        xlim([-100 3000]);
%         ylim([0.6 1.5]);
        ylim([0.3 1.8]);

        title(s.name,'Interpreter', 'none');
        set(gca,'XDir','reverse')

        config_plot(f);

        saveas(f,[ktrFolder,filesep,s.name,'_ktr',savename_suffix,'.png']);
        
        close(f);
    end
    

end
end


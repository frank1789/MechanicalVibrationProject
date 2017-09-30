function graphicssinesweep(pCount, pNamefile, pEstimatedtf, pEstimatedtf2)

if ~exist('pEstimatedtf2','var')|| isempty(pEstimatedtf2)
    figure();
    bodeplot(pEstimatedtf, logspace(0,2,1000));
    grid on;
    namefile = strcat(pNamefile, ['bodediagram' int2str(pCount)]);
%     legend('free damping', 'proportinal damping','Location','west');
%     saveas(gcf,namefile,'epsc')
    
    
    
else   
    figure();
    bodeplot(pEstimatedtf,pEstimatedtf2, logspace(0,2,1000));
    grid on;
    namefile = strcat(pNamefile, ['bodediagram' int2str(pCount)]);
%     legend('free damping', 'proportinal damping','Location','west');
%     saveas(gcf,namefile,'epsc')    
    
end

end
function graphicssinesweep(pCount, pEstimatedtf, pNamefile)

    figure();
    bodeplot(pEstimatedtf, logspace(0,2,1000));
    grid on;
    namefile = strcat(pNamefile, ['bodediagram' int2str(pCount)]);
%     legend('free damping', 'proportinal damping','Location','west');
%     saveas(gcf,namefile,'epsc')
end
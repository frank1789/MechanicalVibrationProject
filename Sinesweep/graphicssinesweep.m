function graphicssinesweep(pCount, pNamefile, pEstimatedtf, pEstimatedtf2)
%graphicssinesweep generate plot of estimate transfer function compute by
%tfest and store images.eps

if ~exist('pEstimatedtf2','var')|| isempty(pEstimatedtf2)
    % print single tf
    figure();
    bodeplot(pEstimatedtf, logspace(0,2,1000));
    grid on;
    str = sprintf('Bode Diagram %s', pNamefile);
    title(str);
    namefile = strcat(pNamefile, ['bodediagram' int2str(pCount)]);
    saveas(gcf,namefile,'epsc');
else
    %print overlap tf
    figure();
    bodeplot(pEstimatedtf,pEstimatedtf2, logspace(0,2,1000));
    grid on;
    namefile = strcat(pNamefile, ['bodediagram' int2str(pCount)]);
    legend('sine sweep slow', 'sine sweep fast');
    saveas(gcf,namefile,'epsc');
end

end
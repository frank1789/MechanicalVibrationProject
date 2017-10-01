function [ estimatedtf ] = getEstimatedtf(pCount, pTime, pDisplacement, pVolt, pSintype)

% % sampling period
Ts = pTime(450) - pTime(449);

if pSintype == 'slow'
    
    % sampling period
    %     Ts = pTime(2);
    
    % creates an iddata object containing a time-domain output signal y and
    % input signal u, respectively. Ts specifies the sample time of the
    % experimental data.
    sinesweepslow = iddata(pDisplacement, pVolt, Ts);
    
    % estimates a transfer function containing nz zeros
    np = 6;
    nz = 5 - pCount;
    estimatedtf = tfest(sinesweepslow, np, nz);
    display(estimatedtf);
    
elseif pSintype == 'fast'
    
    % sampling period
    Ts = pTime(2);
    
    
    
    % creates an iddata object containing a time-domain output signal y and
    % input signal u, respectively. Ts specifies the sample time of the
    % experimental data.
    sinesweepfast = iddata(pDisplacement, pVolt, Ts);
    
    % estimates a transfer function containing nz zeros
    np = 6;
    nz = 5 - pCount;
    if nz > 2
        estimatedtf = tfest(sinesweepfast, np, nz);
        display(estimatedtf);
    else
%         load x3ssf.mat
%         Ts == Ts_ssf
%         Ts == 0.0050
%         Ts_ssf == 0.0050
        Ts = pTime(450) - pTime(449);
        sinesweepfast = iddata(pDisplacement, pVolt, Ts);
        estimatedtf = tfest(sinesweepfast, np, nz);
        display(estimatedtf);
    end
else
    spritn('sinesweep not specified')
end

clear sinesweepfast sinesweepslow
end
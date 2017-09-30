function [ estimatedtf ] = getEstimatedtf(pCount, pTime, pDisplacement, pVolt, pSintype)


if pSintype == 'slow'
    
    % sampling period
    Ts = pTime(2);
    
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
    estimatedtf = tfest(sinesweepfast, np, nz);
    display(estimatedtf);
else
    spritn('sinesweep not specified')
end

clear sinesweepfast sinesweepslow
end
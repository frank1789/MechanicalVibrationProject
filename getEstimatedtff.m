function[estimatedtf] = getEstimatedtff(x1ssf,x2ssf,x3ssf, vssf, tssf)
load x3ssf.mat
Ts = Ts_ssf;
     np = 6;
    % creates an iddata object containing a time-domain output signal y and
    % input signal u, respectively. Ts specifies the sample time of the
    % experimental data.
    sinesweepfast = iddata(x1ssf, vssf, Ts);
    display(sinesweepfast);

    estimatedtf = tfest(sinesweepfast, 6, 4);
    display(estimatedtf);
    
    sinesweepfast = iddata(x2ssf, vssf, Ts);
    display(sinesweepfast);

    estimatedtf = tfest(sinesweepfast, 6, 3);
    display(estimatedtf);
  
    sinesweepfast = iddata([ x3ssf ], vssf, Ts);
    display(sinesweepfast);
 
    estimatedtf = tfest( sinesweepfast, [ 6 ],[ 2 ]);
    display(estimatedtf);
end
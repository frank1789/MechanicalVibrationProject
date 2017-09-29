function [ g_v_per, R31_per, R32_per ] = geterrorspring(pGain_v, pRatio_k3_k2, pRatio_k3_k1, pInputdata)
%geterrorspring compute the ratio between nominal value of spring and ratio
% between nomial gain value and coputed gain value

% compute nomial ratio 
nominalRatio_31 = pInputdata.stiffness.k3 / pInputdata.stiffness.k1;
nominalRatio_32 = pInputdata.stiffness.k3 / pInputdata.stiffness.k2;

% compute error between nominal and estimated spring
R31_per = abs((pRatio_k3_k1 - nominalRatio_31) * 100 / (pRatio_k3_k1));
R32_per = abs((pRatio_k3_k2 - nominalRatio_32) * 100 / (pRatio_k3_k2));
display(R31_per);
display( R32_per);

% compute error between nominal and estimated gain
g_v_per = abs(100* (pInputdata.voltage.v - pGain_v) / (pInputdata.voltage.v));
display(g_v_per);
end
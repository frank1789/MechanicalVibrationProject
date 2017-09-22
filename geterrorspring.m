function [] = geterrorspring(pGain_v, pRatio_k3_k2, pRatio_k3_k1, pInputdata, pDataset)


nominalRatio_31 = pInputdata.stiffness.k3 / pInputdata.stiffness.k1;
nominalRatio_32 = pInputdata.stiffness.k3 / pInputdata.stiffness.k2;


R31_per = abs((pRatio_k3_k1 - nominalRatio_31) * 100 / (pRatio_k3_k1));
R32_per = abs((pRatio_k3_k2 - nominalRatio_32) * 100 / (pRatio_k3_k2));
display(R31_per);display( R32_per);

g_v_per = abs(100* (pGain_v - pDataset.voltage.v(1)) / (pDataset.voltage.v(1)));
display(g_v_per);
end
function [ gain_v, Ratio_k3_k2, Ratio_k3_k1 ] = computeSteadyStateRatioStiff( pK, pDataset, pInputdata)
%computeSteadyStateRatioStiff

% local simbolic value
syms gain_v v k3

% multiply times k_3 supposing is the correct value
gainVoltMatrix = inv(pK) * [gain_v * v; 0; 0];

% multiply the matrix by k3/v, where suppose k3 the best value(correct)
gaVoltMatk3v = expand(simplify(gainVoltMatrix * (k3/v)));
disp('Show matrix for research voltage-to-force:');
display(gaVoltMatk3v);

% instaziate local variable to indetify the interval of time where the responde is
% flat
time = pDataset.time.t(2);

% instaziate local variable to store flat interval at time
x1_coefficient = 0;
x2_coefficient = 0;
x3_coefficient = 0;

% compute flat interval
for i = 0:3
    x1_coefficient = mean(pDataset.Displacement.x1((1.5 / time + i * 10) : (4.5 / time + i * 10)));
    x2_coefficient = mean(pDataset.Displacement.x2((1.5 / time + i * 10) : (4.5 / time + i * 10)));
    x3_coefficient = mean(pDataset.Displacement.x3((1.5 / time + i * 10) : (4.5 / time + i * 10)));
end

% compute the factor stiffness nominal stiffness of spring k3 / nominal voltage
spring_volt = pInputdata.stiffness.k3 / pDataset.voltage.v(1);

% compute from last equation of simbolic "gaVoltMatk3v" the estimated gain
gain_v = x3_coefficient * spring_volt;

% compute the ratio beetwen spring K3/K2, K3/K1
syms ratio_3_2 ratio_3_1

Ratio_k3_k2 = double(solve(gain_v + (gain_v * ratio_3_2) == x2_coefficient * spring_volt, ratio_3_2));
Ratio_k3_k1 = double(solve(gain_v + (gain_v * ratio_3_1) + (gain_v * Ratio_k3_k2) == x1_coefficient * spring_volt, ratio_3_1));
fprintf('The ratio the ratio beetwen spring K3/K2 = %f, K3/K1 = %f\n', Ratio_k3_k2, Ratio_k3_k1);

end
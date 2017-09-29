function [G] = gettranferfunc(pOptimValue, pInputdata)
%gettranferfunc return the transfer funtion for generate bodle plot
%
%How to works:
% instaziate local variable from Optimalvalue and data
% k1 = pInputdata.stiffness.k1;
% k2 = pInputdata.stiffness.k2;
% k3 = pInputdata.stiffness.k3;
% 
% % assemble stifness matrix [K]
% K = [k1 -k1 0; -k1 k2+k1 -k2; 0 -k2 k2+k3];
% 
% % M = eye(3,3);
% % C = zeros(3,3);
% 
% if length(pOptimValue) > 6
%     % assign local variable for mass
%     m1 = pOptimValue(1);
%     m2 = pOptimValue(2);
%     m3 = pOptimValue(3);
%     
%     % assemble mass matrix [M]
%     M = [m1 0 0; 0 m2 0; 0 0 m3];
%     
%     % assign local variable for damping
%     c1 = pOptimValue(4);
%     c2 = pOptimValue(5);
%     c3 = pOptimValue(6);
%     
%     % assemble damping matric [C]
%     C = [c1 0 0; 0 c2 0; 0 0 c3];   
%     
%     % assign local variable for gain scale factor
%     gain = pOptimValue(7);
%     
% else
%     % assign local variable for mass
%     m1 = pOptimValue(1);
%     m2 = pOptimValue(2);
%     m3 = pOptimValue(3);
%     
%     % assemble mass matrix [M]
%     M = [m1 0 0; 0 m2 0; 0 0 m3];
%     
%     % assign local variable for damping
%     alpha = pOptimValue(5);
%     beta = pOptimValue(6);
%     
%     % proportional damping [C] = alpha * [M] + beta * [K]
%     C = alpha * M + beta * K;  
%     
%     % assign local variable for gain scale factor
%     gain = pOptimValue(4);
%     
% end
% 
% s = tf('s');
% A = M * s^2 + C * s + K;
% G = (gain) * tf(inv(A));
% end


% instaziate local variable from Optimalvalue and data
k1 = pInputdata.stiffness.k1;
k2 = pInputdata.stiffness.k2;
k3 = pInputdata.stiffness.k3;

% assemble stifness matrix [K]
K = [k1 -k1 0; -k1 k2+k1 -k2; 0 -k2 k2+k3];

% M = eye(3,3);
% C = zeros(3,3);

if length(pOptimValue) > 6
    % assign local variable for mass
    m1 = pOptimValue(1);
    m2 = pOptimValue(2);
    m3 = pOptimValue(3);
    
    % assemble mass matrix [M]
    M = [m1 0 0; 0 m2 0; 0 0 m3];
    
    % assign local variable for damping
    c1 = pOptimValue(4);
    c2 = pOptimValue(5);
    c3 = pOptimValue(6);
    
    % assemble damping matric [C]
    C = [c1 0 0; 0 c2 0; 0 0 c3];   
    
    % assign local variable for gain scale factor
    gain = pOptimValue(7);
    
else
    % assign local variable for mass
    m1 = pOptimValue(1);
    m2 = pOptimValue(2);
    m3 = pOptimValue(3);
    
    % assemble mass matrix [M]
    M = [m1 0 0; 0 m2 0; 0 0 m3];
    
    % assign local variable for damping
    alpha = pOptimValue(5);
    beta = pOptimValue(6);
    
    % proportional damping [C] = alpha * [M] + beta * [K]
    C = alpha * M + beta * K;  
    
    % assign local variable for gain scale factor
    gain = pOptimValue(4);
    
end

s = tf('s');
A = M * s^2 + C * s + K;
G = (gain) * tf(inv(A));
end
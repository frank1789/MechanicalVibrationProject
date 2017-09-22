function [ F ] = computeforce( pDatatset )
%Computeforce the function compute the force applied at system passed
%argument voltage from dataset and return a matrix with force component
                                                 
% costant from readme
% k_a = 2;                               % servo amp gain
% k_t = 0.1;                             % Servo motor Torque
% k_mp = 26.25;                          % Motor Pinioon pitch radius inverse 

% F = zeros(length(data_steps.time.t), 3);         % initialize F's vector

% temporary vector of force
% f = (Inputdata.k_a * Inputdata.k_t * Inputdata.k_mp) * data_steps.voltage.v;   

% return matrix of force
% F(:,1) = f;                                      

% initialize F's vector
F = zeros(length( pDatatset.time.t), 3);

k_a = 2;                               % servo amp gain
k_t = 0.1;                             % Servo motor Torque
k_mp = 26.25;                          % Motor Pinioon pitch radius inverse 

% tmp vector 
f = (k_a * k_t * k_mp) * pDatatset.voltage.v;

% return matrix of forces
F(:,1) = f;

% free memory
clear f

end
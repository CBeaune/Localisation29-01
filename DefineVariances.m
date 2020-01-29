% Set the parameters which have the "Determined by student" comment to tune
% the Kalman filter. Do not modify anything else in this file.

% Uncertainty on initial position of the robot.

sigmaX     = 0.1;         % Determined by student 
sigmaY     = 0.1;         % Determined by student 
sigmaTheta = 0.1*pi/180 ;   % Determined by student 
Pinit = diag( [sigmaX^2 sigmaY^2 sigmaTheta^2] ) ;


% Measurement noise.

%Constantes nÃ©cessaires au calcul
T = samplingPeriod ; %Periode d echantillonage des capteurs
Vmax = topRobotSpeed ; %Vitesse maximale admissible par le robot.

sigmaMeasurement = (width+T*Vmax)/sqrt(12) ; %determine par nous ; Vmax*T ; %change le 29/01
Qgamma = sigmaMeasurement^2;

% Input noise

sigmaTuning = 0.007 ; %change le 29/01
Qwheels = sigmaTuning^2 * eye(2) ;
Qbeta   = jointToCartesian * Qwheels * jointToCartesian.' ; 

% State noise

Qalpha = zeros(3) ;

% Mahalanobis distance threshold

mahaThreshold = sqrt(chi2inv(0.95,1)); % Determined by student

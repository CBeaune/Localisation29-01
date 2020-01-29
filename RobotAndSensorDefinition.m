%OK

% Constants defining the robot and sensor setup
% All lengths are in mm.

% Robot characteristics

rwheel           = 21.5 ;      % Wheel radius
trackGauge       = 112  ;      % Distance between the fixed wheels
realencoderRes = 180  ;      % In dots per wheel rotation

%Fausser l'odometrie
%Ajouté par Charlotte le 15/01
dumbfactor = 10 ;
encoderRes = realencoderRes/dumbfactor;


% To dumb down the robot's odometry, we artificially lower the
% encoder resolution by a certain factor.
% Also, we use a lower sampling frequency

samplingFrequency    = 100 ;
samplingPeriod       = 1/samplingFrequency ;

% The sensor is supposed to be orthogonal to axis Xm of the robot.

noLineDetected =  0   ;      % Bit value when no magnet is detected
lineDetected   =  1   ;      % Bit value when a magnet is detected

topRobotSpeed = 300 ;%change le 29/01


% Homogeneous coordinates of the line detector sensors in robot frame Rm.

% One column per sensor.
            
%s1   s2  ...

mSensors = [  0   0  ;
             50 -50  ;
              1   1  ] ;

nbSensors = size(mSensors,2) ;                  


% Line spacing   


xSpacing = 305 ;

ySpacing = 305 ;

width = 5 ;


% ---------------------------------------------------------------

% The following are calculated from previous data. Do not modify.

% ---------------------------------------------------------------

hwidth = width/2;

dots2rad = (2*pi)/encoderRes ;

rad2dots = 1/dots2rad        ;


jointToCartesian = [ rwheel/2 rwheel/2 ; rwheel/trackGauge -rwheel/trackGauge ] ;
jointToCartesianfaux = [ (rwheel+0.2)/2 rwheel/2 ; (rwheel+0.2)/trackGauge -rwheel/trackGauge ] ; %change le 29/01

cartesianToJoint = inv(jointToCartesian) ; 

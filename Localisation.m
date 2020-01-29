%main code, a utiliser apres CreateRobotTraj et simulate_sensors

% Localization using a grid of magnets in the ground.
% -----
% Usage: 
%    - Set the characteristics of the robot and sensor in 
%      RobotAndSensorDefinition.m
%    - Set noise levels in DefineVariances.m
%    - Set robot initial position in the present file.
%    - Execute this file.
%    - Then, execute PlotResults.m to get the plots for analysis.
% You can also use, for a brand new fresh start each time:
% clear all; close all; MagnetLoc; PlotResults;
% -----
% Project history:
%    - Project initiator and principle of the sensor: GaÃ«tan Garcia
%    - Sensor manufacturing and test: JoÃ«l Rousseau.
%    - Characterization of the sensor: EMARO1/ARIA1 project, Hendry Chame
% and Marcelo Gaudenzi de Faria. Supervision: GaÃ«tan Garcia
%    - First implementation (Matlab): EMARO1/ARIA1 project, Filip Cichosz
% and Marteen Samuel. Supervision: GaÃ«tan Garcia.
%    - This program (using Samuel and Cichosz's work): GaÃ«tan Garcia

RobotAndSensorDefinition ;
DefineVariances ;

%Load the data file
dataFile = uigetfile('*.mat','Select data file') ;
if isunix 
    load(dataFile);
else
    load(dataFile);
end

X = [ xreal(1), yreal(1), atan2( (yreal(2)-yreal(1)) , (xreal(2)-xreal(1)) ) ].' ;%position intiale
%pour l'instant on prend la vraie position initiale

P = Pinit ; 

% Log data for display of results. Do not modify.

LogData( 0 , 'init' , X , P , [0;0] ) ;

% Skip motionless parts of the data at beginning and end of the experiment
% and return only meaningful data, with wheel rotations in radians.
% Also reduce encoder resolution and frequency according to factors
% set in RobotDefinition.m

wbHandle = waitbar(0,'Computing...') ;

for i = 2 : length(treal) 
    
    t = (i-1)*samplingPeriod ;
    
    waitbar(i/length(treal)) ;

    % Calculate input vector from proprioceptive sensors
    deltaq = [ qR(i) - qR(i-1) ; 
               qL(i) - qL(i-1) ] ;
    U = jointToCartesianfaux * deltaq ;  % joint speed to Cartesian speed.
    
    % Predic state (here odometry)
    X = EvolutionModel( X , U ) ;
    
    % Calculate linear approximation of the system equation
    A = [ 1 0 -U(1)*sin(X(3));
          0 1 U(1)*cos(X(3)) ;
          0  0   1] ;
      
    B = [ cos(X(3)) 0;
          sin(X(3)) 0;
          0 1] ;
       
    % Error propagation
    P = A*P*(A.') + B*Qbeta*(B.') + Qalpha ;
    
    LogData( t , 'prediction' , X , P , U , [0;0] ) ;
    
    % Vector of sensor data at time t
    measures = sensorState(i,:) ; 
            
    %On ignore les mesures proches des angles (on utilise la position
    %reelle) 
    for measNumber = 1 : length(measures) 
        if measures(measNumber)==lineDetected 
            % Homogeneous transform of robot frame with respect to world frame
            oTm = [ cos(X(3)) -sin(X(3)) X(1);
                    sin(X(3))  cos(X(3)) X(2);
                    0   0   1] ;

            % Now in homogeneous coordinates for calculations.
            mMeasSensor = mSensors(:,measNumber);

            % Corresponding position in absolute frame. 
            oMeasSensor = oTm * mMeasSensor ;            
            
            %on determine si x ou y est le plus proche d'un multiple de la
            %taille des carreaux
            oMeasSensorNormalized=oMeasSensor./[xSpacing;ySpacing;1];
            oMeasSensorNormalized=mod(oMeasSensorNormalized,[1,1,1]);
            if oMeasSensorNormalized(1)>0.5
                oMeasSensorNormalized(1)=1-oMeasSensorNormalized(1);
            end
            if oMeasSensorNormalized(2) >0.5
                oMeasSensorNormalized(2)=1-oMeasSensorNormalized(2);
            end
            
            if oMeasSensorNormalized(1)<oMeasSensorNormalized(2) %une verticale est plus proche qu'une horizontale 
                
                % Due to measurement and localization errors, the previously calculated
                % position does not match an vetical line.
                % Which actual vertical line is closest? It will be the candidate line for
                % the update phase of the state estimator.
                oRealSensor=oMeasSensor;
                oRealSensor(1) = round(oMeasSensor(1)/xSpacing)*xSpacing ;

                % The position of the real magnet in robot frame. It will in general 
                % be different from the measured one. 
                mRealSensor = oTm \ oRealSensor ;  % That's inv(oTm)*oRealMagnet = mTo*oRealMagnet

                % Measure: x coordinate of the sensor in Rm.
                Y = mSensors(1,measNumber);

                % The expected measurement is the x coordinate of the real 
                % line in the robot frame.
                Yhat = mRealSensor(1) ;

                %C = [ 1 0 -mSensors(1,measNumber)*sin(X(3))-mSensors(2,measNumber)*cos(X(3))];
                %change pour corriger le filtre
                C = [ 1 0 -X(1)*sin(X(3))-X(2)*cos(X(3))];
                innov = Y - Yhat ;   
                dMaha = innov * sqrt( 1 / ( C*P*C.' + Qgamma) ) ;
                LogData( t , 'measurement' , X , P , [0;0] , mSensors(1:2,measNumber) ) ;

                if dMaha <= mahaThreshold & abs(oMeasSensorNormalized(1)-oMeasSensorNormalized(2))>0.3 %change le 29/01
                    K = 1/( C*P*C.' + Qgamma) * P * C.'  ;
                    X = X + innov.*K ;
                    P = (eye(length(X)) - K*C) * P ;
                    LogData( t , 'update' , X , P , [0;0] , [0;0] ) ;
                end
            else %ligne horizontale
                
               % Due to measurement and localization errors, the previously calculated
                % position does not match an vetical line.
                % Which actual vertical line is closest? It will be the candidate line for
                % the update phase of the state estimator.
                oRealSensor=oMeasSensor;
                oRealSensor(2) = round(oMeasSensor(2)/ySpacing)*ySpacing ;

                % The position of the real magnet in robot frame. It will in general 
                % be different from the measured one. 
                mRealSensor = oTm \ oRealSensor ;  % That's inv(oTm)*oRealMagnet = mTo*oRealMagnet

                % Measure: x coordinate of the sensor in Rm.
                Y = mSensors(2,measNumber);

                % The expected measurement is the x coordinate of the real 
                % line in the robot frame.
                Yhat = mRealSensor() ;

                %C = [ 0 1 mSensors(1,measNumber)*cos(X(3))-mSensors(2,measNumber)*sin(X(3))];
                C = [ 0 1 X(1)*cos(X(3))-X(2)*sin(X(3))];
                innov = Y - Yhat ;   
                dMaha = abs(innov) * sqrt( 1 / ( C*P*C.' + Qgamma) ) ;
                LogData( t , 'measurement' , X , P , [0;0] , mSensors(1:2,measNumber) ) ;

                if dMaha <= mahaThreshold & abs(oMeasSensorNormalized(1)-oMeasSensorNormalized(2))>0.3 %change le 29/01
                    K = 1/( C*P*C.' + Qgamma) * P * C.'  ;
                    X = X + innov.*K ;
                    P = (eye(length(X)) - K*C) * P ;
                    LogData( t , 'update' , X , P , [0;0] , [0;0] ) ;
                end
            end
        end
    end
end

% Save all tunings and robot parameters, so the conditions under
% which the results were obtained are known. Also save inputs and 
% measurements for later display.    
    
save inputLog ...
     rwheel trackGauge encoderRes samplingPeriod ...
     xSpacing ySpacing  ...
     nbSensors U sensorState ...
     Pinit Qgamma sigmaTuning Qbeta Qalpha mahaThreshold 

 
LogData( t , 'termination' , X , P , [0;0] , [0;0] ) ;
close(wbHandle) ;
%close all;
PlotResults;

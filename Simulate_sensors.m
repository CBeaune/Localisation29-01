%OK

RobotAndSensorDefinition ;

%Load the data file
dataFile = uigetfile('*.mat','Select data file') ;
if isunix 
    load(dataFile);
else
    load(dataFile);
end

% Resample at the desired frequency (should be lower than the frequency 
% obtained using the simulation, say ten times less).
% Data to resample: qRight, qleft, xq, yq, theta
totalTime = tq(length(tq)) ;
nbSamples = floor(totalTime/samplingPeriod) ;
treal     = [0:nbSamples]*samplingPeriod ; 
treal     = treal.' ;
qR        = interp1 (tq, qRight, treal) ;
qL        = interp1 (tq, qLeft , treal) ;
xreal     = interp1 (tq, xq    , treal) ;
yreal     = interp1 (tq, yq    , treal) ;
thetareal = interp1 (tq, thetaq, treal) ;

% Apply quantization noise on encoder values based on resolution
qR = round(qR*rad2dots)*dots2rad ;
qL = round(qL*rad2dots)*dots2rad ;

% Now simulate the line detector state along the path, taking into 
% account their location with respect to the robot frame.
% Principle: The ground is assumed to be a checkerboard floor. The sensor
% state is 1 if the sensor is above a white square, 0 otherwise.
% Let xs,ys be the absolute coordinates of the sensor and xSpacing,ySpacing
% the x and y dimensions of the rectangles of the checkerboard floor.
% Assuming the rectangle [0,xSpacing]x[0,ySpacing] is white, then:
% if floor(xs/xSpacing) + floor(ys/xSpacing) is even, the sensor is above 
% a white square, otherwise it's above a black square.

% This is just for checking graphically with two sensors. Not general code.
%nbWhites = [0 0] ;
%nbBlacks = [0 0];

sensorState = zeros( length(treal) , nbSensors ) ;
figure;
for i = 1 : nbSamples
    for j = 1 : nbSensors
        oTm = [ cos(thetareal(i))  ,  -sin(thetareal(i))  ,  xreal(i)  ;
                sin(thetareal(i))  ,   cos(thetareal(i))  ,  yreal(i)  ; 
                      0        ,         0        ,    1   ] ;
        oSensor = oTm * mSensors(:,j) ;
        xs = oSensor(1) ;
        ys = oSensor(2) ;

	%checkboard floor
        % sensorState(i,j) = ~rem( floor(xs/xSpacing)+floor(ys/ySpacing) , 2 ) ;
	%square floor with lines
        sensorState(i,j) = abs( floor(xs/xSpacing)*xSpacing - xs )<hwidth | abs( ceil(xs/xSpacing)*xSpacing - xs )<hwidth | abs( ceil(ys/ySpacing)*ySpacing - ys )<hwidth | abs( ceil(ys/ySpacing)*ySpacing - ys )<hwidth; 
        if sensorState(i,j)
            plot(xs,ys,'g+') ;
            hold on; 
        else
            plot(xs,ys,'r+') ;
            hold on; 
        end
    end
end
for i = 1 : ceil(5000/xSpacing)
    plot(linspace(0,5000,5000),ones(1,5000)*i*xSpacing , 'k--' , 'LineWidth', 0.3 );
    hold on;
    plot(ones(1,5000)*i*xSpacing,linspace(0,5000,5000) , 'k--' , 'LineWidth', 0.3 );
    hold on;
end
save simu dots2rad rad2dots rwheel trackGauge topRobotSpeed ...
          mSensors xSpacing ySpacing ...
          treal xreal yreal thetareal qR qL sensorState ;

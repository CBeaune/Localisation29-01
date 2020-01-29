% A script to display the results of the last execution of MagnetLoc.
% MagnetLoc stores its results in 'log.txt' and the meaningful input data
% in inputLog.mat (among which robot parameters and noise variances).
% Graphs displayed:
% Figure 1: 
%   - The path calculated by odometry only (red).
%   - The path estimated by the Kalman filter (blue).
%   - The real path (green).
%   - The locations of the magnets which have been detected (black dots).
%   - The estimated locations of the detected line in absolute frame
%       using the measurement (that's variable oMeasMagnet in the program).
% Figure 2: 
%   - Speed and rotation speed, as estimated using the encoders.
% Figure 3: 
%   - Estimated error standard deviations (extracted directly from
%       the diagonal of P, hence in absolute frame.
% Figure 4: 
%   - Estimated error standard deviations in robot frame.
% Figure 5: 
%   - Mahalanobis distances calculated with the line closest to 
%       the measurement point (candidate line) in blue.
%   - Mahalanobis distances calculated with the four nearest neighbor 
%       lines of candidate line (2 with the same orientatin and 2 with a different orientation )in red.
%   - Mahalanobis distance threshold used in the program (black line).
% Figure 6: 
%   - Estimated x, y, theta as functions of time.
% Figure 7:
%   - Number of sensor on a line at each time instant.
% Figure 8:
%   - Raw sensor measurements as a function of the curvilinear abscissa
%       (distance traveled by point M). The vertical axis represents the
%       state of each sensor. A vertical line indicates a closed 
%       sensor (a sensor which detects a line).
%   - You may comment out this graph when you don't need it anymore 
%       (when you're done estimating the measurement noise).


% Load the inputs to the problem (robot charateristics + tuning +
% speed and rotation speed + measurements. They have been saved by
% MagnetLoc in inputLog.mat.

if ~exist('inputLog.mat','file')
    disp('File inputLog.mat not found. Did you run MagnetLoc first?');
    return
end
if ~exist('log.txt','file')
    disp('File log.txt not found. Did you run MagnetLoc first?');
    return
end

load inputLog ...
     nbSensors samplingPeriod xSpacing ySpacing  ...
     sensorState Qgamma mahaThreshold

 
% Load the results calculated by MagnetLoc, logged in log.txt.
fid = fopen('log.txt','r');
firstline = fgetl(fid) ;
numvars = numel(strread(firstline,'%s')); % #ok<DSTRRD>
fgetl(fid); %<-- Skip the second line
data = textscan(fid,repmat('%f',1,numvars)); % #ok<NASGU>
% The "unused variable" warning has been suppressed: variable "data" is
% used in the "eval" instruction, but code analyzer does not see it.

% Next instruction sets all variables whose name are on the first line
% of the file, here calcPhase, t, x, y, theta, P11...P33, U1, U2, Y1, Y2.
eval(['[' strrep(firstline,' ',',') '] = deal(data{:});']) ;
fclose(fid);

% Prepare vectors and matrices
nbRes = length(t) ;
nbPeriods = sum(calcPhase==1) ;
U    = zeros(2,nbPeriods) ;
Xodo = zeros(3,nbPeriods+1) ;

% Reconstruct the inputs. Suppress lines that are non prediction, since 
% their values have been set to zero conventionally in the logs.
% Keeping those would be impractical when displaying velocities.
U(1,:) = U1(find(calcPhase==1)) ;
U(2,:) = U2(find(calcPhase==1)) ;
tOdo = t(find(calcPhase==1));

% Calculate the measured positions of the magnets in the absolute frame
% and the magnets closest to these measured positions for display. 
% Do this only when calcPhase = 2 (measurement).
nbMeasurementPhases = sum(calcPhase==2) ;
estLinePos   = zeros(2,nbMeasurementPhases) ;
exactLinePos = zeros(2,nbMeasurementPhases) ;
j = 0 ;
for i = 1 : nbRes 
    if calcPhase(i) == 2
        j = j+1 ;
        oTm = [ cos(theta(i)) , -sin(theta(i)) , x(i) ;
                sin(theta(i)) ,  cos(theta(i)) , y(i) ;
                      0       ,        0       ,  1   ] ;
        oEstimatedMagnet = oTm * [ y1(i) ; y2(i) ; 1 ] ;
        
        oEstimatedSensorNormalized=oEstimatedMagnet./[xSpacing;ySpacing;1];
        oEstimatedSensorNormalized=mod(oEstimatedSensorNormalized,[1,1,1]);
        if oEstimatedSensorNormalized(1)>0.5
            oEstimatedSensorNormalized(1)=1-oEstimatedSensorNormalized(1);
        end
        if oEstimatedSensorNormalized(2) >0.5
            oEstimatedSensorNormalized(2)=1-oEstimatedSensorNormalized(2);
        end
            
         if oEstimatedSensorNormalized(1)<oEstimatedSensorNormalized(2) %une verticale est plus proche qu'une horizontale 
            oExactMagnetPos  = round( oEstimatedMagnet ./ [xSpacing ; 1 ; 1] ) .* [xSpacing ; 1 ; 1] ;
         else 
            oExactMagnetPos  = round( oEstimatedMagnet ./ [1 ; ySpacing ; 1] ) .* [1 ; ySpacing ; 1] ;
         end
          
        estLinePos(:,j)   = oEstimatedMagnet(1:2) ;
        exactLinePos(:,j) = oExactMagnetPos(1:2)  ;
    end
end

% Compute odometry only estimated path
Xodo(:,1) = [x(1) ; y(1) ; theta(1)] ;
for i = 1 : nbPeriods
    Xodo(:,i+1) = Xodo(:,i) + ...
        [ U(1,i)*cos(Xodo(3,i)) ;
          U(1,i)*sin(Xodo(3,i)) ;
          U(2,i) ] ;
end
travDistance = zeros(1,nbPeriods) ;
for i = 2 : nbPeriods
    travDistance(i) = travDistance(i-1) + U(1,i) ;
end

% Plot robot path, Kalman filter estimation and odometry only estimation

figure; 
plot( Xodo(1,:), Xodo(2,:) , 'r' , 'LineWidth', 2 ) ;
hold on ;
plot( x,y , 'b' , 'LineWidth', 2 ) ;
hold on ;
plot( xreal,yreal , 'k' , 'LineWidth', 1 ) ;
zoom on ; grid on; axis('equal');
title('Real path (black), estimated path EKF (blue) and odometry (red)');
xlabel('x (mm)');
ylabel('y (mm)');

% On top of the path, indicate estimated and real magnet positions.

hold on;
plot( estLinePos(1,:), estLinePos(2,:) , 'gx' ) ;
hold on;
plot( exactLinePos(1,:), exactLinePos(2,:) , 'kx' ) ;
hold on
for i = 1 : length(exactLinePos(1,:))
    if mod(exactLinePos(1,i),xSpacing)==0
        plot( ones(1,5000)*exactLinePos(1,i),linspace(0,5000,5000) , 'k--' , 'LineWidth', 0.3 );
    else
        plot(linspace(0,5000,5000),ones(1,5000)*exactLinePos(2,i) , 'k--' , 'LineWidth', 0.3 );
    end
end

% Plot odometry-estimated speed and rotation speed

figure; 
subplot(2,1,1);
plot( tOdo,U(1,:)/samplingPeriod , 'LineWidth',2 );
xlabel('t (s)');
ylabel('v (mm/s)');
title('Odometry-estimated speed');
zoom on ; grid on;
subplot(2,1,2);
plot( tOdo,U(2,:)*180/pi/samplingPeriod , 'LineWidth',2 );
xlabel('t (s)')
ylabel('w (deg/s)' , 'LineWidth',2 );
title('Odometry-estimated rotation speed');
zoom on ; grid on;


% Plot estimated variances in absolute reference frame

sigx     = sqrt(P11) ;
sigy     = sqrt(P22) ;
sigtheta = sqrt(P33) ;
tRes = t( find(calcPhase==3) ) ;
figure;
subplot(3,1,1);
maximum = max(sigx) ;
for k=1:numel(tRes) 
    line([tRes(k) tRes(k)],[0 maximum],'Color','g','LineStyle',':' ,...
       'LineWidth',2 );      
end
hold on ;
plot( t,sigx , 'LineWidth',2 );
xlabel('t (s)') ;
ylabel('sigma_x (mm)') ;
title('Estimated standard deviations in absolute ref. frame');
title('Est. std dev. in abs. ref. frame');
zoom on ; grid on;
set(gca,'FontSize',14)
subplot(3,1,2);
maximum = max(sigy) ;
for k=1:numel(tRes) 
   line([tRes(k) tRes(k)],[0 maximum],'Color','g','LineStyle',':' , ...
       'LineWidth',2 );      
end
hold on ;
plot( t,sigy , 'LineWidth',2 );
xlabel('t (s)') ;
ylabel('sigma_y (mm)');
zoom on ; grid on;
subplot(3,1,3);
maximum = max(sigtheta*180/pi) ;
for k=1:numel(tRes) 
   line([tRes(k) tRes(k)],[0 maximum],'Color','g','LineStyle',':' ,...
       'LineWidth',2 );      
end
hold on ;
plot( t,sigtheta*180/pi , 'LineWidth',2 );
xlabel('t (s)') ;
ylabel('sigma_{theta} (deg.)');
zoom on ; grid on;

% Calculate covariance matrix in frame Rm
msigx     = zeros(1,length(t)) ;
msigy     = zeros(1,length(t)) ;
msigtheta = zeros(1,length(t)) ;
for i = 1 : length(t)
    m_Omega_o = [  cos(theta(i))  , -sin(theta(i))  ,  0  ; %changement de signe des sin 
                   sin(theta(i))  ,  cos(theta(i))  ,  0  ;
                        0         ,       0         ,  1  ] ;
    oP = [ P11(i)  ,  P12(i)  ,  P13(i)  ;
           P12(i)  ,  P22(i)  ,  P23(i)  ;
           P13(i)  ,  P23(i)  ,  P33(i)  ] ;
    mP = m_Omega_o * oP * m_Omega_o.' ;
    msigx(i)     = sqrt( mP(1,1) )   ;
    msigy(i)     = sqrt( mP(2,2) )   ;
    msigtheta(i) = sqrt( mP(3,3) )   ;
end

% Plot variances in robot frame.

figure;
subplot(3,1,1);
maximum = max(msigx) ;
for k=1:numel(tRes) 
   line([tRes(k) tRes(k)],[0 maximum],'Color','g','LineStyle',':' ,...
       'LineWidth',2 );      
end
hold on ;
plot( t,msigx , 'LineWidth',2 );
xlabel('t (s)') ;
ylabel('sigma_x (mm)');
title('Estimated standard deviations in robot frame');
zoom on ; grid on;
subplot(3,1,2);
maximum = max(msigy) ;
for k=1:numel(tRes) 
   line([tRes(k) tRes(k)],[0 maximum],'Color','g','LineStyle',':' , ...
       'LineWidth',2 );      
end
hold on ;
plot( t,msigy , 'LineWidth',2 );
xlabel('t (s)') ;
ylabel('sigma_y (mm)');
zoom on ; grid on;
subplot(3,1,3);
maximum = max(msigtheta*180/pi) ;
for k=1:numel(tRes) 
   line([tRes(k) tRes(k)],[0 maximum],'Color','g','LineStyle',':' ,...
       'LineWidth',2 );      
end
hold on ;
plot( t,msigtheta*180/pi , 'LineWidth',2 );
xlabel('t (s)') ;
ylabel('sigma_{theta} (deg.)');
zoom on ; grid on;

% Calculate Mahalanobis distances, including for lines that are the 
% closest neighbors of the line closest to measurement point.   

tLineDetection = zeros(1,sum(calcPhase==2)) ;
dMahaAll = zeros(5,sum(calcPhase==2)) ;
j = 0 ;
%change le calcul de C dans le rep�re absolu %change le 29/01
for i = 1 : length(t) 
    
    if calcPhase(i) ~= 2 
        continue ;  % Not a measurement phase
    end
     
    j = j+1 ;
    tLineDetection(j) = t(i) ;
    
    % Calculate homogeneous transform of the robot with respect to the world frame
    oTm = [ cos(theta(i)) , -sin(theta(i)) , x(i)  ;
            sin(theta(i)) ,  cos(theta(i)) , y(i)  ;
                  0       ,        0       ,  1    ] ;
    
    % Measurement vector in homogeneous coordinates for calculations.  The
    % measurement vector hould be a scalar
    mMeasSensor = [ y1(i) ; y2(i) ; 1 ] ;
    
    % Corresponding position in absolute frame
    oMeasSensor = oTm * mMeasSensor ;
    
    % Which actual line is closest to the estimated position?
    oMeasSensorNormalized=oMeasSensor./[xSpacing;ySpacing;1];
    oMeasSensorNormalized=mod(oMeasSensorNormalized,[1,1,1]);
    if oMeasSensorNormalized(1)>0.5
        oMeasSensorNormalized(1)=1-oMeasSensorNormalized(1);
    end
    if oMeasSensorNormalized(2) >0.5
        oMeasSensorNormalized(2)=1-oMeasSensorNormalized(2);
    end

     if oMeasSensorNormalized(1)<oMeasSensorNormalized(2) %une verticale est plus proche qu'une horizontale 
        oRealSensor  = round( oMeasSensor ./ [xSpacing ; 1 ; 1] ) .* [xSpacing ; 1 ; 1] ;
            
        % The position of the real magnet in robot frame
        mRealMagnet = oTm \ oRealSensor ; % That's inv(oTm)*oRealMagnet
        
        %Measuerment vector
        Y=mMeasSensor(1:2); %we need both coordinates because even if candidate line is vertical, neighbor line could be vertical or horizontal
        % The expected measurement are the two coordinates of the real
        % magnet in the robot frame.
        Yhat = mRealSensor(1) ;
        C = [ 1 0 -oMeasSensor(1)*sin(theta(i))-oMeasSensor(2)*cos(theta(i))];

        innov = Y(1) - Yhat ;
        P = [ P11(i)  ,  P12(i)  ,  P13(i)  ;
              P12(i)  ,  P22(i)  ,  P23(i)  ;
              P13(i)  ,  P23(i)  ,  P33(i)  ] ;
        dMaha = abs(innov) * sqrt( 1 / ( C*P*C.' + Qgamma) ); 

        dMahaAll(1,j)       = dMaha            ;
        estLinePos(:,j)   = oMeasSensor(1:2) ;
        exactLinePos(:,j) = oRealSensor(1:2) ;

        % Offset vectors to generate the neighbors, in homogeneous coordinates.

        deltas = [ xSpacing  -xSpacing      0           0       ;
                      0          0       ceil(y(i)/ySpacing)*ySpacing-y(i)   -floor(y(i)/ySpacing)*ySpacing+y(i)   ;
                      0          0          0           0       ] ;

        for neighborIndex = 1 : 4
            oPneighbor = oRealSensor + deltas(:,neighborIndex) ;

            % The position of the line in robot frame is the expected measurement
            if neighborIndex<=2 % vertical neighbors
                YhatNeighbor = oTm \ oPneighbor ; % That's inv(oTm)*oPneighbor
                CNeighbor = [ 1 0 -oMeasSensor(1)*sin(theta(i))-oMeasSensor(2)*cos(theta(i))];
                innovNeighbor = Y(1) - YhatNeighbor(1) ;    % Not in homogeneous coordinates.
                dMahaNeighbor = abs(innovNeighbor) * sqrt( 1 / ( CNeighbor*P*CNeighbor.' + Qgamma) ) ;
                dMahaAll(neighborIndex+1,j) = dMahaNeighbor ;
                
            else %horizontal neighbors
                YhatNeighbor = oTm \ oPneighbor ; % That's inv(oTm)*oPneighbor
                CNeighbor = [ 0 1 oMeasSensor(1)*cos(theta(i))-oMeasSensor(2)*sin(theta(i))];
                innovNeighbor = Y(2) - YhatNeighbor(2) ;    % Not in homogeneous coordinates.
                dMahaNeighbor = abs(innovNeighbor) * sqrt( 1 / ( CNeighbor*P*CNeighbor.' + Qgamma) ) ;
                dMahaAll(neighborIndex+1,j) = dMahaNeighbor ;
            end
        end
        
     else %closest is horizontal
        oRealSensor  = round( oMeasSensor ./ [ySpacing ; 1 ; 1] ) .* [ySpacing ; 1 ; 1] ;
            
        % The position of the real magnet in robot frame
        mRealMagnet = oTm \ oRealSensor ; % That's inv(oTm)*oRealMagnet
        
        %Measuerment vector
        Y=mMeasSensor(1:2); %we need both coordinates because even if candidate line is vertical, neighbor line could be vertical or horizontal
        % The expected measurement are the two coordinates of the real
        % magnet in the robot frame.
        Yhat = mRealSensor(2) ;
        C = [ 0 1 oMeasSensor(1)*cos(theta(i))-oMeasSensor(2)*sin(theta(i))];

        innov = Y(2) - Yhat ;
        P = [ P11(i)  ,  P12(i)  ,  P13(i)  ;
              P12(i)  ,  P22(i)  ,  P23(i)  ;
              P13(i)  ,  P23(i)  ,  P33(i)  ] ;
        dMaha = abs(innov) * sqrt( 1 / ( C*P*C.' + Qgamma) ) ;

        dMahaAll(1,j)       = dMaha            ;
        estLinePos(:,j)   = oMeasSensor(1:2) ;
        exactLinePos(:,j) = oRealSensor(1:2) ;

        % Offset vectors to generate the neighbors, in homogeneous coordinates.

        deltas = [ ceil(x(i)/xSpacing)*xSpacing-x(i)  -floor(x(i)/xSpacing)*xSpacing+x(i)      0           0       ;
                      0          0       ySpacing   -ySpacing   ;
                      0          0          0           0       ] ;

        for neighborIndex = 1 : 4
            oPneighbor = oRealSensor + deltas(:,neighborIndex) ;

            % The position of the line in robot frame is the expected measurement
            if neighborIndex<=2 % vertical neighbors
                YhatNeighbor = oTm \ oPneighbor ; % That's inv(oTm)*oPneighbor
                CNeighbor = [ 1 0 -oMeasSensor(1)*sin(theta(i))-oMeasSensor(2)*cos(theta(i))];
                innovNeighbor = Y(1) - YhatNeighbor(1) ;    % Not in homogeneous coordinates.
                dMahaNeighbor = abs(innovNeighbor) * sqrt( 1 / ( CNeighbor*P*CNeighbor.' + Qgamma) ) ;
                dMahaAll(neighborIndex+1,j) = dMahaNeighbor ;
                
            else %horizontal neighbors
                YhatNeighbor = oTm \ oPneighbor ; % That's inv(oTm)*oPneighbor
                CNeighbor = [ 0 1 oMeasSensor(1)*cos(theta(i))-oMeasSensor(2)*sin(theta(i))];
                innovNeighbor = Y(2) - YhatNeighbor(2) ;    % Not in homogeneous coordinates.
                dMahaNeighbor = abs(innovNeighbor) * sqrt( 1 / ( CNeighbor*P*CNeighbor.' + Qgamma) ) ;
                dMahaAll(neighborIndex+1,j) = dMahaNeighbor ;
            end     
        end
     end
    
end

% Plot Mahalanobis distances. Blue dots are for closest magnet,
% red dots are for neighbor magnets.

figure; 
plot( tLineDetection , dMahaAll(1,:) , 'bo' , 'LineWidth',2 ) ;
for k = 2:3
    hold on; 
    plot( tLineDetection , dMahaAll(k,:) , 'ro' , 'LineWidth',2 ) ;
end
for k = 4:5
    hold on; 
    plot( tLineDetection , dMahaAll(k,:) , 'yo' , 'LineWidth',2 ) ;
end
hold on;
plot( tLineDetection , mahaThreshold*ones(1,size(dMahaAll,2)) , 'k' ,...
    'LineWidth',2 ) ; 
xlabel('t (s)');
ylabel('Mahalanobis distance (no dimension).');
title({'Mahalanobis distances:','closest line (blue), vertical neighbors (red) and horizontal neighbors (yellow)'});
zoom on; grid on; 

% Plot x, y and theta as functions of time

figure; 

subplot(3,1,1);
plot( t,x , 'LineWidth',2 );
xlabel('t (s)')
ylabel('x (mm)');
title('Position and heading as functions of time.');
zoom on ; grid on;

subplot(3,1,2);
plot( t,y , 'LineWidth',2 );
xlabel('t (s)')
ylabel('y (mm)');
zoom on ; grid on;

subplot(3,1,3);
plot( t,theta*180/pi , 'LineWidth',2 );
xlabel('t (s)')
ylabel('theta (deg.)');
zoom on ; grid on;

% Determine the number of measurements (i.e. the number of sensor on a line)
% at each step.

for time=1:length(treal)
    Mesures=0;
    for sensor = 1:nbSensors
        Mesures=Mesures+sensorState(time,sensor);
    end
    nbMeas(time)=Mesures;
end

% Plot number of measurements at each time step.

figure; 
plot(treal,nbMeas,'o' , 'LineWidth', 2 ) ;
xlabel('time (s)');
%yticks([0 1 2]);
title('Number of lines detected at each step.');
zoom on; grid on;


% Plot raw sensor measurements
rawMeas = zeros( nbSensors , nbPeriods ) ;
for k = 1 : nbPeriods
    rawMeas( : , k ) = sensorState(k,:) ;
end
figure; 
for n = 1 : nbSensors
    for k = 1 : nbPeriods
        if rawMeas(n,k) == 1
            hold on ;
            line([travDistance(k) travDistance(k)],[n-0.5 n+0.5],'Color','b','LineStyle','-' , 'LineWidth', 2 );
        end
    end
end
set(gca,'YLim',[0 nbSensors+1]) ;
xlabel('Travelled distance of point M (mm)');
ylabel('Sensor number') ;
title('State of sensors (blue = line dectected)') ;
zoom on; grid on;

% Plot errors in absolute frame.
%Fait par Sandro le 15/01

figure;
subplot(2,1,1);
xref=zeros(size(t));
yref=zeros(size(t));
thetaref=zeros(size(t));
c=0;
troissigma = zeros(size(t)); %ajouté par Audrey
troissigmay = zeros(size(t)); %ajouté par Audrey
for i=1:length(t)
    if calcPhase(i)==0 | calcPhase(i)==1
        xref(i)=xreal(i-c);
        yref(i)=yreal(i-c);
        thetaref(i)=thetareal(i-c);
    else
        xref(i)=xref(i-1);
        yref(i)=yref(i-1);
        thetaref(i)=thetaref(i-1);
        c=c+1;
    end
    troissigma(i)=3*sigmaMeasurement; %ajouté par Audrey
end
plot( t,troissigma,'k', 'LineWidth',1 );
hold on;
plot( t,-troissigma,'k', 'LineWidth',1 );
hold on;
plot( t,x-xref ,'b', 'LineWidth',2 );
hold on;
xlabel('t (s)') ;
ylabel('erreur_x (mm)');
title('Error in absolute frame');
zoom on ; grid on;

subplot(2,1,2);
plot( t,troissigma,'k', 'LineWidth',1 );
hold on;
plot( t,-troissigma,'k', 'LineWidth',1 );
hold on;
plot( t,y-yref ,'b', 'LineWidth',2 );
hold on;
xlabel('t (s)') ;
ylabel('erreur_y (mm)');
zoom on ; grid on;

% subplot(3,1,3);
% plot( t,3*sigtheta*180/pi ,'k', 'LineWidth',1 );
% hold on ;
% plot( t,-3*sigtheta*180/pi ,'k', 'LineWidth',1 );
% hold on ;
% plot( t,theta-thetaref ,'b', 'LineWidth',2 );
% hold on;
% xlabel('t (s)') ;
% ylabel('sigma_{theta} (deg.)');
% zoom on ; grid on;


% Calculate and display odometry error.

fprintf('\nTotal travelled distance: %d mm\n',round(sum(abs(U1))));
fprintf('Final odometry error: %3.1f %%\n\n', ...
    (norm([xreal(size(xreal,1));yreal(size(yreal,1))]-Xodo(1:2,size(Xodo,2))) / sum(abs(U1)) )*100 );

% Calculate and display KF error.

fprintf('\nTotal travelled distance: %d mm\n',round(sum(abs(U1))));
fprintf('Final KF error: %3.1f %%\n\n', ...
    (norm([xreal(size(xreal,1));yreal(size(yreal,1))]-X(1:2,size(X,2))) / sum(abs(U1)) )*100 );

% Calculate percentage of rejected closest lines:

fprintf('Lines rejected: %3.1f %%\n', ...
    100*numel(find(dMahaAll(1,:) > mahaThreshold ))/numel(dMahaAll(1,:)));

fprintf('Neighbor lines under threshold: %3.1f %%\n\n', ...
    100*numel(find(dMahaAll(2:5,:) <= mahaThreshold ))/numel(dMahaAll(2:5,:))) ;

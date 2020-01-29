%OK

%%Create a trajectory with the help of splines that connect points defined by the user

RobotAndSensorDefinition ;

axis([0 5000 0 5000])
for i = 1 : floor(5000/xSpacing)
    plot( [i*xSpacing i*xSpacing] , [0 5000] ) ;
    hold on;
    plot( [0 5000] , [i*xSpacing i*xSpacing] ) ;
    hold on;
end
    
% Initially, the list of points is empty.
x = [] ;
y = [] ;
n = 0;
% Loop, picking up the points.
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
while but == 1
    [xi,yi,but] = ginput(1) ;
    plot(xi,yi,'ro')
    x = [ x ; xi ] ;
    y = [ y ; yi ] ;
end

% define waypoints
t = [0:length(x)-1];
% calculate spline for way points
tq = 0:0.001:length(x)-1;
tq = tq.' ;
slope0 = 0;
slopeF = 0;
xq = spline(t,[slope0; x; slopeF],tq);
yq = spline(t,[slope0; y; slopeF],tq);

% plot spline in t-x, t-y and x-y space

% Plot the waypoints
figure; 
plot(x,y,'ro') ;

% Plot the interpolated curve.
hold on; 
plot(xq,yq,'b.');
hold off;

topSpeed = 0 ;
for i=2 : length(tq)
    v(i-1) = norm([xq(i);yq(i)]-[xq(i-1);yq(i-1)])/(tq(i)-tq(i-1)) ;
    topSpeed = max(topSpeed,v(i-1)) ;
end

% Scale time so the top speed has a given value.
speedFactor = topRobotSpeed/topSpeed;
tq = tq/speedFactor ;

% Calculate robot orientation at each point. Done numericall, but 
% can be done using the spline coefficients. I just don't bother...
thetaq(1) = atan2( yq(2)-yq(1) , xq(2)-xq(1) ) ;
for i=2 : length(tq)-1
    thetaq(i) = atan2( yq(i+1)-yq(i) , xq(i+1)-xq(i) ) ;
    % Ensure continuity of the orientation
    while thetaq(i)-thetaq(i-1) > pi 
        thetaq(i) = thetaq(i)-2*pi ;
    end
    while thetaq(i)-thetaq(i-1) < -pi 
        thetaq(i) = thetaq(i)+2*pi ;
    end 
end
thetaq(length(tq)) = thetaq(length(tq)-1) ;
thetaq = thetaq.' ;

% Calculate elementary travelled distance and elementary rotation of the
% robot between time instants.
deltaD     = zeros(length(tq),1) ;
deltaTheta = zeros(length(tq),1) ;
for i=2 : length(tq)
    deltaD(i) = norm([xq(i);yq(i)]-[xq(i-1);yq(i-1)]) ;
    deltaTheta(i) = thetaq(i)-thetaq(i-1) ;
end

% Calculate wheel rotation angle at each time instant (zero at t=0)
qRight = zeros(length(tq),1) ;
qLeft  = zeros(length(tq),1) ;
for i = 2 : length(tq) 
    delta = cartesianToJoint*[deltaD(i);deltaTheta(i)] ;
    qRight(i) = qRight(i-1) + delta(1) ;
    qLeft (i) = qLeft (i-1) + delta(2) ;  
end

% Check results by recalculating path using odometry equations.
xOdo     = zeros(length(tq),1) ;
yOdo     = xOdo ; 
thetaOdo = xOdo ;

xOdo(1)     = xq(1)     ;
yOdo(1)     = yq(1)     ;
thetaOdo(1) = thetaq(1) ;
for i = 2 : length(tq) 
    dCart = jointToCartesian*[ qRight(i)-qRight(i-1) ; qLeft(i)-qLeft(i-1) ] ;
    xOdo(i)     = xOdo(i-1)     + dCart(1)*cos(thetaOdo(i-1)) ;
    yOdo(i)     = yOdo(i-1)     + dCart(1)*sin(thetaOdo(i-1)) ;
    thetaOdo(i) = thetaOdo(i-1) + dCart(2)                    ;
end 

figure ; 
plot(xq,yq,'b','LineWidth',4) ;
hold on; 
plot(xOdo,yOdo,'r','LineWidth',2) ;
axis('equal') ;

% Save all robot information and the trajectory information
save traj rwheel trackGauge tq xq yq thetaq qRight qLeft

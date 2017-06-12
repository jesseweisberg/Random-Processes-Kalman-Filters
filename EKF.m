clear all;
%%%
% Simulated Extended Kalman Filter for ballistic vehicle tracking
% Author: Jesse Weisberg
% 
% vehicle dynamical model: y'' = (r)*(y')^2 - g, r= constant air density factor
%
% written as state-space model:
% x1 = y, x2 = y' (states: position, velocity)
% x' = f(x) + Gu = f(x) + [0 Su]'u
% 
% measurements of vertical positions are taken directly (discrete
% intervals):
% z_k = y(t_k)+v_k, E[v_k*v_k'] = R_k (measurement variance)
% 
% f, F, phi(tau), phi_k, H were determined by hand, 
%%%

% constants   
Su = 1.5;    %Adjustable noise term for Q (tuning the EKF)
             % Su=1.5 just about minimizes rms errors (mean errors ~=0,
             % position variance ~=36, which is slightly larger than actual
             % rms errors, velocity errors behaved well too)
Sp = 1000;   %Std. dev. of initial position
Sv = 100;    %Std. dev. of initial velocity
Sm = 100;    %Std. dev. of measurement error
y0 = 100000; %Initial average height
v0 = -500;   %Initial average velocity
 
vt = 1500;    %terminal velocity
g  = 32.2;    %acceleration due to gravity
r  = g/vt^2;  %rho
 
R   = Sm^2;  %measurement variance
T   = 1;     %sample period
Max = 100;   %Number of steps in simulation
 
InitHeight = random('Normal',y0,Sp);  %True initial height
InitVel    = random('Normal',v0,Sv);   %True initial velocity
h          = -atanh(InitVel/vt);
 
%Series to be plotted. (plot errors and filter variance terms.) 
Time  = 0*(1:Max);
TruthPos = 0*(1:Max);
TruthVel = 0*(1:Max);
EstimatedPos = 0*(1:Max);
EstimatedVel = 0*(1:Max);
SPosError = 0*(1:Max);
SuError = 0*(1:Max);
Simulation = 0*(1:Max);
 
%System matrices 
I   = eye(2);
H = [1 0]; 
P = [Sp^2   0   ;   %initial P
       0   Sv^2]; 
X   = [y0 ;
       v0];
 
%Each loop here processes a measurement.   
for k = 1:Max
   t       = k*T;                                      %Time
   Time(k) = t;
   v = - vt*tanh(g*t/vt + h);                          %True velocity
   y = y0 - vt^2/g*(log(cosh(g*t/vt+h))-log(cosh(h))); %True height
   Z = y + random('Normal',0,Sm);                      %Measured height
   TruthPos(k) = y; 
   TruthVel(k) = v;
        
   K = P*H'/(H*P*H'+R); %Kalman gain
   X = X + K*(Z - H*X); %State correction
   P = (I-K*H)*P;
   
   %Phi, Q matrices
   Phi = [1 T;
          0 1+2*T*r*v];
      
   Q = [(1/3)*(Su^2*T^3) (2/3)*(T^3*Su^2*r*v) + (1/2)*(Su^2*T^2);
        (2/3)*(T^3*Su^2*r*v) + (1/2)*(Su^2*T^2) Su^2*((4/3)*r^2*v^2*T^3 + 2*r*v*T^2+T)];
    
   EstimatedPos(k) = X(1);%Covariance correction
   EstimatedVel(k) = X(2);  
   
   % collect performance statistics
   SPosError(k) = P(1,1)^.5;  

   f = [v;
        r*v^2-g];
   X = X + f*T;         %State prediction
   P = Phi*P*Phi' + Q;  %Covariance prediction
   %Simulation(k) = X(1);
   
end
%At this point, processed 100 seconds of measurements.
sumPos = 0;
sumVel = 0;
sumVelS = 0;
sumPosS = 0;

for k = 50:Max
    VelError = TruthVel(k) - EstimatedVel(k);
    PosError = TruthPos(k) - EstimatedPos(k);
    VelErrorS = (TruthVel(k) - EstimatedVel(k))^2;
    PosErrorS = (TruthPos(k) - EstimatedPos(k))^2;
    sumVel = sumVel + VelError;
    sumPos = sumPos + PosError;  
    sumVelS = sumVelS + VelErrorS;
    sumPosS = sumPosS + PosErrorS;
end
meanVelError = sumVel/50;
meanPosError = sumPos/50;

RMSVelError = (sumVelS/50)^.5;
RMSPosError = (sumPosS/50)^.5;

disp(meanVelError);
disp(meanPosError);
disp(RMSVelError);
disp(RMSPosError);

% plot position error and velocity error
plot(Time, SPosError, Time, -SPosError, Time, TruthPos-EstimatedPos) 
% plot(Time, SPosError, Time, -SPosError, Time, TruthVel-EstimatedVel)

%could calculate cumulative statistics (e.g., rms) here.



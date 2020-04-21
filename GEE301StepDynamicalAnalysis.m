timeinterval = 0.01;
time = 0:timeinterval:2;


%% MASS SPRING DAMPER PROPERTIES
m = 8.5;
ctest = 300;
ktest = 13000;
bodymass = 60;

%%MASS SPRING DAMPER STATE SPACE MODELLING
% State space explanation - https://www.youtube.com/watch?v=hpeKrMG-WP0
% state space matlab function - http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=SystemModeling
MSDA = [0 1; -ktest/m -ctest/m];
MSDB = [0 1/m]';
MSDC = [1 0];
MSDD = [0];

MSDsys = ss(MSDA,MSDB,MSDC,MSDD);
[n,d] = ss2tf(MSDA,MSDB,MSDC,MSDD);
MSDtf = tf(n,d);

%%INPUT FORCE AT EVERY TIME INTERVAL
%trapezoidal waveform underneath -https://www.google.com/search?q=ground+reaction+force+stairs&rlz=1C1CHBF_en-GBGB822GB822&tbm=isch&source=iu&ictx=1&fir=a4e9jGXnVEXmqM%253A%252CPVG7hBcPRbTVbM%252C_&usg=AI4_-kQjBZUucZlSe8-wYw8IiDjnQ4ydEg&sa=X&ved=2ahUKEwjIucqMz7TgAhWysHEKHX_FA8sQ9QEwAXoECAUQBA#imgrc=a4e9jGXnVEXmqM:
f(1,1:((2/timeinterval)+1)) = 0;
f1(1,1:((2/timeinterval)+1)) = 0;
%initial peak
for i = 1:((0.1/timeinterval)+1)
 f(1,i+1) = -(9.81*(1.4*bodymass/10)*i);
end
%force decrease in step
for i = 1:19
 f(1,11+i) = -(9.81*1.4*bodymass - 9.81*(1.737)*i);
end
%Constant for section
for i = 1:30
 f(1,30+i) = -(9.81*bodymass*0.85);
end
%Rapid fallaway to 0
for i = 1:10
 f(1,60+i) = -(9.81*bodymass*0.85 - 50.031*i);
end

%initial condition of displacement and velocity
x0 = [-0.05 0];

%built-in dynamical analysis tool
y = lsim(MSDsys,f,time);

%velocity function
v(1,1:((2/timeinterval)+1)) = 0;
for i = 1:(length(y)-1)
    distancedifference = y(i+1,1) - y(i,1);
    timedifference = timeinterval;
    v(1,i+1) = distancedifference/timedifference;
end

%% DC GENERATOR PROPERTIES - http://www.moog.com/literature/MCG/moc23series.pdf
%% C34-L70-10
% J = 0.0003; %moment of inertia of rotor
% b = 0.05; %motor viscous friction (peak torque/rated rotational speed)
% Ke = 0.041; %electromotive force constant (rated torque/rated current)
% Kt = Ke;%motor torque constant
% R = 0.14; %electrical resistance of armature
% L = 0.24; %electrical inductance of coil
% r = 0.05; %gear radius
% 
% %%DC GENERATOR STATE SPACE MODELLING
% %http://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=SystemModeling
% %http://www2.nkfust.edu.tw/~tuky/servo/PDF/Ch2.pdf
% DCGA = [-R/L -Ke/L; Ke/J -b/J];
% DCGB = [1/L 0; 0 -1/J];
% DCGC = [1 0];
% DCGD = [0 0];
% 
% voltageapplied(1,1:((2/timeinterval)+1)) = 0;
% torqueapplied(1,1:((2/timeinterval)+1)) = 0;
% 
% for i = 1:((2/timeinterval)+1)
%     torqueapplied(1,i) = ctest*v(1,i)*r;
% end
% 
% DCGsys = ss(DCGA,DCGB,DCGC,DCGD);
% % DCGsystf=tf(DCGsys);
% % DCGVoltageInputTF = DCGsystf(1,1);
% % DCGTorqueInputTF = DCGsystf(1,2);
% DCGinput = [voltageapplied; torqueapplied];
% CurrentOutDCG = lsim(DCGsys, DCGinput, time); %output is current
% 
% TorqueOutDCG = CurrentOutDCG.*Ke;

subplot(322)
plot(time,y);
grid on
title('Displacement of step over time');

subplot(321)
plot(time,f);
grid on
title('Force Regime');

subplot(323)
plot(time,v);
grid on
title('Velocity of Step over time')

subplot(324)
vsq = v.*v;
power = vsq.*ctest;
plot(time, power);
title('Power extracted by the damper');
energy = trapz(time,power);


DamperForce = power./(-v);
subplot(325)
plot(time,DamperForce);
title('Force exerted by the Damper');
% 
% subplot(326)
% plot(time,TorqueOutDCG);


%% check energy in damper is correct using damper value and velocity
% Check method found using first pdf in https://www.google.com/search?q=energy+dissipated+in+damper&rlz=1C1CHBF_en-GBGB822GB822&oq=energy+dissipated+in+damper&aqs=chrome..69i57j0.10625j1j4&sourceid=chrome&ie=UTF-8
% 'Module: 8 Lecture: 1 Energy dissipated by damping - nptel' - downloaded
% on laptop
vsq = v.*v;
CHECKenergyfactor = trapz(time,vsq);
CHECKenergy = CHECKenergyfactor*ctest;

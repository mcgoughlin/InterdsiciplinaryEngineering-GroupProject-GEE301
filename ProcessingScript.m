clear

load GearRatio2ohm4_1InertiaAdjusted2Data
time = 0:0.01:2;
PowerRaw = PowerOut.data;
BackTorqueRaw = ElectroTorque.data;
RotorSpeedRaw = RotorSpeed.data;

timeinterval = 0.01;
time = 0:timeinterval:2;
maxdisp = -0.025;

SampleRateT = round(length(BackTorqueRaw)/length(time));
SampleRateP = round(length(PowerRaw)/length(time));
SampleRateR = round(length(RotorSpeedRaw)/length(time));

for i = 1:length(time)-1
    BackTorque(i) = BackTorqueRaw(SampleRateT*i);
    
end
BackTorque(length(time)) = 0;

for i = 1:length(time)-1
    Power(i) = PowerRaw(SampleRateP*i);
end
Power(length(time)) = 0;

for i = 1:length(time)-1
    RotorS(i) = RotorSpeedRaw(SampleRateR*i);
    
end
RotorS(length(time)) = 0;


Energy = trapz(time,Power);
clear PowerOut.mat
clear PowerRaw
clear BackTorqueRaw
%Gear system properties
maxturnsratio = 4.1;
req1 = 0.0205;

plot(PowerOut)
% MASS SPRING DAMPER PROPERTIES
m = 140;
ctest =50;
ktest = 1000;
bodymass = 60;

%%MASS SPRING DAMPER STATE SPACE MODELLING
% State space explanation - https://www.youtube.com/watch?v=hpeKrMG-WP0
% state space matlab function - http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=SystemModeling
MSDA = [0 1; -ktest/m -ctest/m];
MSDB = [0 1/m]';
MSDC = [1 0];
MSDD = [0];

MSDsys = ss(MSDA,MSDB,MSDC,MSDD);

%%INPUT FORCE AT EVERY TIME INTERVAL
%trapezoidal waveform underneath -https://www.google.com/search?q=ground+reaction+force+stairs&rlz=1C1CHBF_en-GBGB822GB822&tbm=isch&source=iu&ictx=1&fir=a4e9jGXnVEXmqM%253A%252CPVG7hBcPRbTVbM%252C_&usg=AI4_-kQjBZUucZlSe8-wYw8IiDjnQ4ydEg&sa=X&ved=2ahUKEwjIucqMz7TgAhWysHEKHX_FA8sQ9QEwAXoECAUQBA#imgrc=a4e9jGXnVEXmqM:
f(1,1:((2/timeinterval)+1)) = 0;
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
x0 = [0 0];

%Finds Input torque to put into Generator Model

inputtorque = f.*(req1/maxturnsratio);

for i = 1:length(f)
    %Turning Back Torque into a Back Force
    inputforce(i) = f(i) - BackTorque(i)*(maxturnsratio/req1);
end

%built-in dynamical analysis tool
y = lsim(MSDsys,inputforce,time);

plot(time,y)

for i = 1:length(y)
    if y(i) < maxdisp
        y(i) = maxdisp;
    end
end

%% This section computes the response of a system that disengages the drivetrain 
%  due to the sprocket freewheel.
m = 10.7;
ctest1 = 20;
ktest1 = ktest;
%different damping value (much lower) as it no longer engages generator


%%MASS SPRING DAMPER STATE SPACE MODELLING
% State space explanation - https://www.youtube.com/watch?v=hpeKrMG-WP0
% state space matlab function - http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=SystemModeling
MSDA1 = [0 1; -ktest1/m -ctest1/m];
MSDB1 = [0 1/m]';
MSDC1 = [1 0];
MSDD1 = [0];

MSDsys1 = ss(MSDA1,MSDB1,MSDC1,MSDD1);

%%INPUT FORCE AT EVERY TIME INTERVAL
%trapezoidal waveform underneath -https://www.google.com/search?q=ground+reaction+force+stairs&rlz=1C1CHBF_en-GBGB822GB822&tbm=isch&source=iu&ictx=1&fir=a4e9jGXnVEXmqM%253A%252CPVG7hBcPRbTVbM%252C_&usg=AI4_-kQjBZUucZlSe8-wYw8IiDjnQ4ydEg&sa=X&ved=2ahUKEwjIucqMz7TgAhWysHEKHX_FA8sQ9QEwAXoECAUQBA#imgrc=a4e9jGXnVEXmqM:

%finds first point where returning spring force = applied force
%This is more accurate when using a finer time interval
springforce = ktest1*maxdisp;
for i = 1:length(f)
    if f(i) == springforce && y(i) == maxdisp
        returntime = time(i);
        break
    elseif f(i) < springforce && f(i+1) > springforce && y(i) == maxdisp
        returntime = time(i+1);
        break
    else
        continue
    end
end

%Applies the 0 discontinuity to the input torque (for the generator model)
%as no more torque can be supplied to the generator when max displacement occurs.
for i = 1:length(inputtorque)
   if y(i) <= maxdisp
       inputtorque(i) = 0;
   else 
       continue
   end
end
% returntime = 0.69;
%initiates matrix for the recovering step motion
f1size = length(time)-returntime/timeinterval;
f1 =zeros(1,f1size);
time1 = zeros(1,f1size);

%ensures force regime is still acting as the step attempts to recover
for i = 1:length(f1)
    f1(i) = f(i+returntime/timeinterval);
end
for i = 1:length(f1)
    time1(i) = time(i+returntime/timeinterval);
end

%initial condition of displacement and velocity
x0 = [y(returntime/timeinterval) 0];

%built-in dynamical analysis tool
y1 = lsim(MSDsys1,f1,time1,x0);

plot(time1,y1);

for i = 1:length(y1)
    if y1(i) < maxdisp
        y1(i) = maxdisp;
    end
end

for i = 1:length(y1)
    if y1(i) > 0
        ReimpactPoint = i;
        for j = ReimpactPoint:length(y1)
            y1(j) = 0;
        end
        break
    end
end


for i = 1:f1size
    y(returntime/timeinterval+i) = y1(i);
end


%velocity function

for i = 1:(length(y)-1)
    distancedifference = y(i+1,1) - y(i,1);
    timedifference = timeinterval;
    v(1,i+1) = distancedifference/timedifference;
end

%find time of max displacement
for i = 1:length(y)
    if y(i) == maxdisp
        maxdisptime = i*timeinterval;
        break
    end
end

%find piezoelectrics reaction forces in terms of time indexes of 0.01s

for i = 1:length(inputforce)
    if i < maxdisptime/timeinterval
        PiezoForce(i) = 0;
    elseif i > returntime/timeinterval
        PiezoForce(i) = 0;
    else
        PiezoForce(i) = f(i);
    end
end


subplot(222)
plot(time,y);
grid on
title('Displacement of step over time');
ylabel(' Displacement / m');
xlabel('Time / s');

subplot(221)
plot(time, f);
grid on
title('Force Regime');
ylabel(' Force / N');
xlabel('Time / s');


subplot(223)
plot(time,v); hold on
plot(time,((RotorS/maxturnsratio)*-req1)); hold off
grid on
title('Velocity of Step over time')
ylabel(' Velocity / m/s');
xlabel('Time / s');

subplot(224)
plot(PowerOut)
ylabel('Power / W')
xlabel('Time / s');

energy = trapz(time,Power);


clear bodymass
clear CHECKenergyfactor
clear ctest
clear ctest1
clear distancedifference
clear f1
clear i
clear f1size
clear j
clear ktest
clear ktest1
clear m
clear maxdisp
clear springforce
clear time1
clear timedifference
clear y1
clear ReimpactPoint
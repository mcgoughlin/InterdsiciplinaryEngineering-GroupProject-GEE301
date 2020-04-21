clear

load ValidationExpData
expdata = expdata';
expdata(431) = 0;
timeinterval = 0.005;
time = 0:timeinterval:2.15;

shaftradius = 0.015;
j = 0.000556;
timeconstant = 0.74;
c = j/timeconstant;

% MASS SPRING DAMPER PROPERTIES
m = j/(shaftradius^2);
btest = c/(shaftradius^2);
ktest = 0.00;

%%MASS SPRING DAMPER STATE SPACE MODELLING

MSDA = [0 1; -ktest/m -btest/m];
MSDB = [0 1/m]';
MSDC = [0 1];
MSDD = [0];

MSDsys = ss(MSDA,MSDB,MSDC,MSDD);

%%INPUT FORCE AT EVERY TIME INTERVAL
%trapezoidal waveform underneath -https://www.google.com/search?q=ground+reaction+force+stairs&rlz=1C1CHBF_en-GBGB822GB822&tbm=isch&source=iu&ictx=1&fir=a4e9jGXnVEXmqM%253A%252CPVG7hBcPRbTVbM%252C_&usg=AI4_-kQjBZUucZlSe8-wYw8IiDjnQ4ydEg&sa=X&ved=2ahUKEwjIucqMz7TgAhWysHEKHX_FA8sQ9QEwAXoECAUQBA#imgrc=a4e9jGXnVEXmqM:
%initial peak
f = zeros(1,length(time));
for i = 9:123
 f(i) = 1.25*9.81;
end


%initial condition of displacement and velocity
x0 = [0 0];

%built-in dynamical analysis tool
y = lsim(MSDsys,f,time);
% for i = 1:(length(y)-1)
%     distancedifference = y(i+1,1) - y(i,1);
%     timedifference = timeinterval;
%     v(1,i+1) = distancedifference/timedifference;
% end

wrads = y./shaftradius;
wrpm = wrads.*(60/(2*3.1416));

figure
plot(time, wrpm);
hold on
plot(time,expdata.*42);
hold off
xlabel('Time / s')
ylabel('Rotational Speed / rpm' )
legend('Equivalent Mass-Spring-Damper System' , 'Experimental Results')


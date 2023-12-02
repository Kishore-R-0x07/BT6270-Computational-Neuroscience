%THIS PROGRAM DEMONSTRATES HODGKIN HUXLEY MODEL IN CURRENT CLAMP EXPERIMENTS AND SHOWS ACTION POTENTIAL PROPAGATION
%Time is in secs, voltage in mvs, conductances in m mho/mm^2, capacitance in uF/mm^2

% threshold value of current is 0.0223

%ImpCur=input('enter the value of the impulse current in microamperes: ');
%TimeTot=input('enter the time for which stimulus is applied in milliseconds');

i_cur = 1;
i_step = 0.01;
imax = 0.6;
firing_rate = zeros(1, int32(imax/i_step));
for ImpCur = 0:i_step:imax

gkmax=.36;
vk=-77; 
gnamax=1.20;
vna=50; 
gl=0.003;
vl=-54.387; 
cm=.01; 

dt=0.01;
niter=50000;
t=(1:niter)*dt;
iapp=ImpCur*ones(1,niter);
%for i=1:100
 %   iapp(1,i)=ImpCur;
 %end;
v=-64.9964;
m=0.0530;
h=0.5960;
n=0.3177;

gnahist=zeros(1,niter);
gkhist=zeros(1,niter);
vhist=zeros(1,niter);
mhist=zeros(1,niter);
hhist=zeros(1,niter);
nhist=zeros(1,niter);


for iter=1:niter
  gna=gnamax*m^3*h; 
  gk=gkmax*n^4; 
  gtot=gna+gk+gl;
  vinf = ((gna*vna+gk*vk+gl*vl)+ iapp(iter))/gtot;
  tauv = cm/gtot;
  v=vinf+(v-vinf)*exp(-dt/tauv);
  alpham = 0.1*(v+40)/(1-exp(-(v+40)/10));
  betam = 4*exp(-0.0556*(v+65));
  alphan = 0.01*(v+55)/(1-exp(-(v+55)/10));
  betan = 0.125*exp(-(v+65)/80);
  alphah = 0.07*exp(-0.05*(v+65));
  betah = 1/(1+exp(-0.1*(v+35)));
  taum = 1/(alpham+betam);
  tauh = 1/(alphah+betah);
  taun = 1/(alphan+betan);
  minf = alpham*taum;
  hinf = alphah*tauh;
  ninf = alphan*taun;
  m=minf+(m-minf)*exp(-dt/taum);
  h=hinf+(h-hinf)*exp(-dt/tauh);
  n=ninf+(n-ninf)*exp(-dt/taun);
  vhist(iter)=v; mhist(iter)=m; hhist(iter)=h; nhist(iter)=n;
end

%figure(1)
%subplot(2,1,1)
%plot(t,vhist)
%title('voltage vs time')

%figure(2)
%subplot(2,1,2)
%plot(t,mhist,'y-', t,hhist,'g.',t,nhist,'b-')
%legend('m','h','n')

%figure(3)
gna=gnamax*(mhist.^3).*hhist; 
gk=gkmax*nhist.^4;
%clf
%plot(t,gna,'r');
%hold on
%plot(t,gk,'b');
%legend('gna','gk')
%hold off

% Finding peaks that cross a threshold value
% And classifying them as Action Potential peaks
% Storing the number of peaks at every ImpCur
thresh = 10;
firing_rate(i_cur) = detect_peaks(thresh, vhist);
i_cur = i_cur + 1;
end

iext = 0:i_step:imax;
% disp(firing_rate);
% Plotting values of I1, I2, I3 too
xline(0.023,'r');
hold on
xline(0.06,'r');
xline(0.45,'r');
plot(iext, firing_rate)
title('Spiking frequency vs Input current')
xlabel('Input current (in uA)')
ylabel('Frequency of Spiking')
hold off

function num_peaks = detect_peaks(thresh, data)
num_peaks = 0;
for p = drange(2:length(data)-1)
    if (data(p) >= thresh) && (data(p) > data(p-1)) && (data(p) > data(p+1))
        num_peaks = num_peaks + 1;
    end
end
end
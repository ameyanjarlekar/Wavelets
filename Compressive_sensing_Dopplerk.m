% this code demonstrates how the doppler spectrum can be estimated using
% CS theory

% the radar return echo, after matched filtering and range-estimation,
% which is used as input to the doppler-estimation block, can be modelled
% as a complex exponential with frequency=doppler frequency

% here, only single PRF is considered for simplicity

% let us now take the CS to Tx side and check the effect

clearvars;
close all;
clc;
%% initialisation %%
nTgt=2;                                 % #targets

Ak=rand(nTgt,1);                                    % signal amplitude
Npulse=5000;                              % #pulses
%prf=500;                                % pulse repeatation frequency
nvar=0.03;                                 % noise variance
m=500;
Rmax=700e3;                             % maximum unambiguous range
prf=5000;                                % pulse repeatation frequency
fdk=round(prf/10*abs(rand(nTgt,1)-0.5));       % target doppler                                 % target doppler
fmk=round(prf*abs(rand(nTgt,1)-0.5));       % target doppler                                 % target doppler
con=round(1*abs(rand(nTgt,1)));       % target doppler                                 % target doppler

%% formulate return echo %%
% selection matrix
A = eye(Npulse);
k = sort(randi(Npulse,[m,1]),'ascend');
A = A(k,:);
ka = (1:Npulse)';
% slow time instants
prt=1/prf;
tm = k*prt;
tma= ka*prt;
PHI1 = zeros(Npulse,nTgt);

% pure received signal
% PHI1=exp(1j*2*pi*tm*fmk');
for i = 1:nTgt
PHI1(:,i)=exp(1j*(2*pi*tma*fmk(i,1)+cos(2*pi*tma*fdk(i,1))));
end
r = PHI1*Ak;
%% noisy measurements %%
w=wgn(Npulse,1,nvar,'linear');
y = r+w;
ra = zeros(Npulse+10,1);
ra(1+5:Npulse+5,1) = y;
fin = y.*conj(ra(1:Npulse));
send = fin(k);
%% sparse recovery %%
% A=2*double(rand(m,Npulse)<=0.5)-1;   % Rademacher random matrix
% A = randn(m,Npulse);                             % gaussian random matrix
% A = orth(A')';
IDFTmat=(1/sqrt(Npulse))*exp(1j*2*pi*(0:Npulse-1)'*(0:Npulse-1)/Npulse);
A1 = A*IDFTmat;
% A1 = myorth(A1.').';
% w=wgn(m,1,nvar,'linear');
% y = A*r + w;

[xhat,S] = OMP(send,A1,0.03,500);

%% plot %%
Fsize=20;
Lwidth=2;
PHIa=exp(1j*cos(2*pi*tma*fdk'));
ra =PHIa*Ak;
figure(1),
set(figure(1),'Name','fft_obtained_from_compressed_Sensing', 'paperpositionmode','auto','paperorientation','landscape');
plot(abs((xhat/max(xhat))))
raf = fft(ra);
figure(2),
set(figure(2),'Name','fft_of_signal_if_only_micro_doppler_is_considered', 'paperpositionmode','auto','paperorientation','landscape');
plot(abs(raf/max(raf)))

ray = (fft(fin));
figure(3),
set(figure(3),'Name','fft_of_correlated_signal', 'paperpositionmode','auto','paperorientation','landscape');
plot(abs(ray/max(ray)))
% fig1=figure(1);
% set(fig1,'Name','fig1', 'paperpositionmode','auto','paperorientation','landscape');
% fk = linspace(-prf/2,prf/2,Npulse);
% plot(fk,abs(fftshift(xhat/max(xhat))),fdk,Ak,'p','Linewidth',Lwidth); grid on;
% xlabel('frequency','Fontsize',Fsize); ylabel('magnitude','Fontsize',Fsize);
% legend('Estimated','Actual')
% set(gca,'Fontsize',Fsize)

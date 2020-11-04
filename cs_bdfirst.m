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
tic
count2 = 0;
count1 = 0;
max_iter = 100;
for iter = 1:max_iter
nTgt=1;                                 % #targets
Ak=rand(nTgt,1);                                    % signal amplitude
while min(Ak,[],'all') < 0.5
    Ak=rand(nTgt,1);     
end
Npulse=1000;                              % #pulses
%prf=500;                                % pulse repeatation frequency
nvar=0.1;                                 % noise variance
m=200;
Rmax=700e3;                             % maximum unambiguous range
prf=1000;                                % pulse repeatation frequency


fmk=round(prf/20*abs(rand(nTgt,1)));        % target doppler                                 % target doppler
if nTgt >1
set = 1;
while set == 1   
   fmk=round(prf/20*abs(rand(nTgt,1))); 
   set = 0;
    for i = 1:nTgt
        for j = 1:nTgt
            if i ~= j
                if (abs(fmk(i,1)-fmk(j,1)) < 10) | (min(fmk,[],'all') < 15)
                    set = 1;
                end
            end
        end
    end
end
else
    while (min(fmk,[],'all') < 15)
        fmk=round(prf/20*abs(rand(nTgt,1))); 
    end
end


fdk=round(prf/4*abs(rand(nTgt,1)));       % target doppler                                 % target doppler
set = 1;
while set == 1   
   fdk=round(prf/4*abs(rand(nTgt,1)));
   set = 0;
    for i = 1:nTgt
        for j = 1:nTgt
            if i ~= j
                if abs(fdk(i,1)-fdk(j,1)) < 100 | (min(fdk,[],'all') < 50)
                    set = 1;
                end
            end
        end
    end
end
            
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
% PHI1=exp(1j*2*pi*tm*fdk');
for i = 1:nTgt
%PHI1(:,i)=exp(1j*(2*pi*tma*fdk(i,1)+cos(2*pi*tma*fmk(i,1))));
%PHI1(:,i)=cos(2*pi*tma*fdk(i,1)+cos(2*pi*tma*fmk(i,1)))+1.0*cos(2*pi*tma*fdk(i,1));
PHI1(:,i)=exp(1j*(2*pi*tma*fdk(i,1)+1.0*cos(2*pi*tma*fmk(i,1))))+1.0*exp(1j*(2*pi*tma*fdk(i,1)));
end
r = PHI1*Ak;
%% noisy measurements %%
snr = 10;
y = awgn(r,snr,'measured');
send = y(k);
%%
%w=wgn(Npulse,1,nvar,'linear');
%y = r+w;
%%%% sparse recovery %%
% A=2*double(rand(m,Npulse)<=0.5)-1;   % Rademacher random matrix
% A = randn(m,Npulse);                             % gaussian random matrix
% A = orth(A')';
IDFTmat=(1/sqrt(Npulse))*exp(1j*2*pi*(0:Npulse-1)'*(0:Npulse-1)/Npulse);
A1 = A*IDFTmat;
% A1 = myorth(A1.').';
% w=wgn(m,1,nvar,'linear');
% y = A*r + w;

[xhat,S] = OMP(send,A1,0.03,500);
% figure(1),
% % set(fig1,'Name','fft_obtained_from_compressed_Sensing', 'paperpositionmode','auto','paperorientation','landscape');
% plot(fftshift(abs((xhat/max(xhat)))))
% comp = abs(fft(y));
% figure(2),
% % set(fig2,'Name','fft_obtained_from_compressed_Sensing', 'paperpositionmode','auto','paperorientation','landscape');
% plot(fftshift(abs((comp/max(comp)))))
%% recon
% figure(3),
% plot(abs(y))
y_recon = IDFTmat*xhat;
% figure(4),
% plot(abs(IDFTmat*xhat))
% figure(5),
% plot(abs(IDFTmat*xhat)-abs(y))
% figure(1),
% plot(fftshift(abs(xhat)))
xhat = fftshift(abs(xhat));
[B,It] = maxk(xhat,nTgt*10);
% I = zeros(2*nTgt,1);
% count1 = 0;
% count2 = 0;
% j = 1;
% for i=1:nTgt*10
%     if j > 4
%         break;
%     end
%     if It(i)>500 & count2<nTgt
%         count2 = count2+1;
%         I(j,1) = It(i);
%     end
%     if It(i)<500 & count1<nTgt
%         count1 = count1+1;
%         I(j,1) = It(i); 
%     end
%     j = j+1;
% end    
% I = sort(I);
ind = zeros(nTgt,1);
ind = abs(It(1:nTgt)-500);
% for i = 1:nTgt
%     ind(i) = ((500-I(i))+(I(2*nTgt+1-i)-500))/2;
% end
fdk = sort(fdk);
ind = sort(ind);
ind = ind-1;
if max(abs(fdk-ind),[],'all') < 5
    count1 = count1 + 1;
end
%%
% ind = sort(ind);
% if nTgt == 2
%     ind (2,1) = ind(2,1)-ind(1,1);
% end
% if nTgt == 3
% ind(3,1) = ind(3,1)-ind(2,1);
% ind (2,1) = ind(2,1)-ind(1,1);
% ind (3,1) = ind(3,1)-ind(1,1);
% end
shift = 7;
fmk_est = zeros(nTgt,1);
for z=1:nTgt
ra = zeros(Npulse+5,1);
ra(1:Npulse,1) = exp(1j*(2*pi*tma*ind(z,1)));
%fin2 = y_recon.*conj(ra(1:Npulse));
fina = y_recon.*conj(ra(1:Npulse));
rb = zeros(Npulse+5+shift,1);
rb(1+shift:Npulse+shift,1) = fina;
fin2 = fina.*conj(rb(1:Npulse));
detect2 = fftshift(abs(fft(fina)));
% figure(2),
% plot(detect2)
p = 250;
[B,I] = maxk(detect2,p);
% I(I>560) = 0;
% I(I<440) = 0;
flag = 0;
for i=1:p
    if flag == 1
        break;
    end
    if I(i) > 0 & abs(500-I(i)) > 10
        for j=1:p
            if flag == 1
                break;
            end
            if abs( abs(I(i)-500)-abs(I(j)-500)) < 4 & abs(I(j)-500) < 55 & (min((abs(abs(500-I(j))-fmk_est(:,1))),[],'all') > 8)
                for l = 1:p
                    if abs(abs(I(l))-2*abs(I(i))) < 5 
                        fmk_est(z,1) = abs(I(i)-500);
                        flag = 1;
                        break;
                    end
                end
            end
            end
        end
    end    
end
fmk = sort(fmk);
fmk_est=fmk_est+1;
fmk_est = sort(fmk_est);
if max(abs(fmk-fmk_est),[],'all') < 3
    count2 = count2 + 1;
end
end
fin2 = count2/max_iter
fin1 = count1/max_iter

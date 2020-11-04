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
count = 0;
max_iter = 100;
for iter = 1:max_iter
nTgt=1;                                 % #targets
Ak=rand(nTgt,1);                                    % signal amplitude
while min(Ak,[],'all') < 0.5
    Ak=rand(nTgt,1);     
end
Npulse=1000;                              % #pulses
%prf=500;                                % pulse repeatation frequency
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
%PHI1(:,i)=cos(2*pi*tma*fdk(i,1)+cos(2*pi*tma*fmk(i,1)));
PHI1(:,i)=exp(1j*(2*pi*tma*fdk(i,1)+cos(2*pi*tma*fmk(i,1))))+1.0*exp(1j*(2*pi*tma*fdk(i,1)));
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
comp = abs(fft(y));
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
%% shift1
% shift = 10;
% ra = zeros(Npulse+shift+5,1);
% ra(1+shift:Npulse+shift,1) = y;
% fin1 = y.*conj(ra(1:Npulse));
%% shift2
shift = 7;
ra = zeros(Npulse+shift+5,1);
ra(1+shift:Npulse+shift,1) = y_recon;
%fin2 = y_recon.*conj(ra(1:Npulse));
fin2 = y_recon.*ra(1:Npulse);

% detect1 = fftshift(abs(fft(fin1)));
Fsize=20;
Lwidth=2;
% figure(1),
% %set(figure(1),'Name','fft_obtained_from_compressed_Sensing', 'paperpositionmode','auto','paperorientation','landscape');
% plot(detect1/max(detect1))
PHIa=cos(cos(2*pi*tma*fmk'));
ra =PHIa*Ak;
raf = fft(ra);
% figure(2),
% %set(figure(2),'Name','fft_of_signal_if_only_micro_doppler_is_considered', 'paperpositionmode','auto','paperorientation','landscape');
% plot(abs(raf/max(raf)))

detect2 = fftshift(abs(fft(fin2)));
Fsize=20;
Lwidth=2;
% figure(3),
% %set(figure(3),'Name','fft_obtained_from_compressed_Sensing', 'paperpositionmode','auto','paperorientation','landscape');
% plot(detect2/max(detect2))
% fig1=figure(1);
% set(fig1,'Name','fig1', 'paperpositionmode','auto','paperorientation','landscape');
% fk = linspace(-prf/2,prf/2,Npulse);
% plot(fk,abs(fftshift(xhat/max(xhat))),fmk,Ak,'p','Linewidth',Lwidth); grid on;
% xlabel('frequency','Fontsize',Fsize); ylabel('magnitude','Fontsize',Fsize);
% legend('Estimated','Actual')
% set(gca,'Fontsize',Fsize)
% %%
% arr = abs(ray/max(abs(ray)));
% avg = mean(arr)
% arr(arr<2*avg) = 0;
% figure(4),
% set(figure(4),'Name','fft_of_correlated_signal', 'paperpositionmode','auto','paperorientation','landscape');
% plot(fftshift(arr))

[B,I] = maxk(detect2,200);
arr_copy = zeros(1000,1);
arr_copy(I) = detect2(I);
% figure(4),
% %set(figure(4),'Name','fft_of_correlated_signal', 'paperpositionmode','auto','paperorientation','landscape');
% plot((arr_copy))

[p q]= size(I);
I = 500-I;
fm_est = ones(nTgt,1)*0;
flag =0;
for i=1:p
    if flag == nTgt
        break
    end
    k = abs(I(i));
    for j = 1:p
        if j ~= i
            if (abs(abs(I(j))-k) < 3) & (abs(I(j)) < 55)  & (min((abs(abs(I(i))-fm_est(:,1))),[],'all') > 8)
                for z = 1:p
                    if (abs(I(z)-2*k) < 3) | (abs(I(z)-3*k) < 3)
                    flag = flag+1;
                    fm_est(flag,1) = abs(I(i));     
                    j = p;
                    break
                    end
                end
            end
        end
    end
end
fm_est;

fm_est = sort(fm_est);
fmk = sort(fmk);
if max(abs(fm_est-fmk),[],'all') < 4
    count = count + 1;
end
end
fin = count/max_iter
toc
% if max(abs(fm_est-fmk),[],'all') < 6
%     count = count + 1;
% end
% end
% count
%%
% [p q]= size(I);
% fm_est = 500;
% flag =0;
% for i=1:p
%     if flag == 1
%         break
%     end
%     k = I(i);
%     l = 1000-k;
%     for j = 1:p
%         if j ~= i
%             if abs(I(j)-l) < 3
%                 flag = flag+1;
%                 fm_est = I(i);
%                 break
%             end
%         end
%     end
% end
% fm_est = abs(500-fm_est);
% %%
% conv1 = exp(-1j*(cos(2*pi*tma*fm_est)));
% out1 = conv1.*fin2;
% findsecond = fftshift(abs(fft(out1)));
% [B2,I2] = maxk(findsecond,20);
% arr_copy = zeros(1000,1);
% arr_copy(I) = detect2(I2);
% figure(5),
% set(figure(5),'Name','fft_of_correlated_signal', 'paperpositionmode','auto','paperorientation','landscape');
% plot((arr_copy))
% [p q]= size(I2);
% fm_est2 = 500;
% flag =0;
% for i=1:p
%     if flag == 1
%         break
%     end
%     k = I2(i);
%     l = 1000-k;
%     for j = 1:p
%         if j ~= i
%             if abs(I2(j)-l) < 3
%                 flag = flag+1;
%                 fm_est2 = I(i);
%                 break
%             end
%         end
%     end
% end
% fm_est2 = abs(500-fm_est2);

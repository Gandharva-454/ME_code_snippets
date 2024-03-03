clc;
clear all;
close;
N = 128; %N=32;
ntap=10;
G = 32;  %%G=8;
nn=N/G;
p=10;
M=4;   %%% QAM constellation
input = randi([0 1], p*100*G, 1);

cp=N/4;
modu1=comm.RectangularQAMModulator(M,'BitInput',true,'NormalizationMethod','Average power', 'AveragePower',1);
modu2=comm.RectangularQAMModulator(M,'BitInput',true,'NormalizationMethod','Average power', 'AveragePower',2.5,'PhaseOffset',(4*pi/16));
modu3=comm.RectangularQAMDemodulator(M,'BitOutput',true,'NormalizationMethod','Average power', 'AveragePower',1);
modu4=comm.RectangularQAMDemodulator(M,'BitOutput',true,'NormalizationMethod','Average power', 'AveragePower',2.5,'PhaseOffset',(4*pi/16));
h=(1/sqrt(ntap*2))*(randn(1,ntap)+sqrt(-1)*randn(1,ntap)); %impulse response of channel
H = fft(h,N); %frequency response of channel
LT=[1 2 3 4;2 3 1 4;3 4 1 2;1 4 2 3];
nt = step(modu1,[0 0 0 1 1 1 1 0]');
lt = step(modu2,[0 0 0 1 1 1 1 0]');
constellation_low = nt;
constellation_high = lt;
modulated_input1=step(modu1,input);
modulated_input2=step(modu2,input);
% modulated_input3=step(modu3,input);
% modulated_input4=step(modu4,input);
scatterplot(modulated_input1);
grid on;
scatterplot(modulated_input2);
grid on;
% scatterplot(modulated_input3);
% grid on;
% scatterplot(modulated_input4);
% grid on;

MASTER = 10*G;   %%% should be a multiple of G
c=0;
for snr=0:1:25
    c = c+1;
    no_symbols = MASTER;
    input = randi([0 1], p*no_symbols*G, 1);
    t = reshape(input,p*G,no_symbols);
    d = zeros(N,no_symbols);
    for i=1:no_symbols
        x=reshape(t(1:(p*G),i),p,G);
        for k=1:G
            p2=x(1:2,k);
            index_dec = bi2de(p2.');
            low_power = LT(index_dec+1,1:2);
            high_power = LT(index_dec+1,3:4);
            subcarr_low = step(modu1, x(3:6,k));
            subcarr_high = step(modu2, x(7:10,k));
            x1 = zeros(nn,1);
            x1(low_power) = subcarr_low;
            x1(high_power) = subcarr_high;
            d((k-1)*nn+1:nn*k,i) = x1;
        end
    end
    ifft_in = reshape(d(:),N,no_symbols);
    input_transform=sqrt(N)*ifft(ifft_in,N,1);
    tx=[input_transform(N-cp+1:N,:)' input_transform']';
    y=[];
    clear i k;
    for i=1:no_symbols
        y=[y.' filter(h,1,tx(:,i)).'].';
    end
    
    s = awgn(y,snr,'measured');
% %     s = y;
    clear i k kk;
    %%%s = y;
    r1=reshape(s,N+cp,no_symbols);
    clear s;
    r2=r1(cp+1:end,:);
    
    r5=sqrt(1/N)*fft(r2,N,1);
    
    r3 = reshape(r5(:),N,no_symbols);
    H1 = repmat(H.',1, no_symbols);
    
    symbol_rec=zeros(1,nn);
    vector_add =[];
    for ii=1:no_symbols
        for k=1:G
            clear cc ll  ERR b;
            cc = r3((k-1)*nn+1:k*nn,ii);
            for kk=1:nn
                clear ERR b;
                for ll = 1:4
                    ERR(ll) = abs(cc(kk) - H1((k-1)*nn+kk,ii)*constellation_low(ll)).^2;
                end
                [b c_low((k-1)*nn+kk)] = min(ERR);
                clear ERR b;
                for ll = 1:4
                    ERR(ll) = abs(cc(kk) - H1((k-1)*nn+kk,ii)*constellation_high(ll)).^2;
                end
                [b c_high((k-1)*nn+kk)] = min(ERR);
            end
        end
        
        clear k;
        for k=1:G
            clear cc int_c_low int_c_high Hnew ERRRR;
            cc = r3((k-1)*nn+1:k*nn,ii);
            int_c_low = c_low((k-1)*nn+1: k*nn);
            int_c_high = c_high((k-1)*nn+1 : k*nn);
             Hnew = H1((k-1)*nn+1 : k*nn,ii);
            %%%%   LT decides the indices of Hnew
            ERRRR(1) = abs(cc(1) - Hnew(1)*constellation_low(int_c_low(1))).^2 + abs(cc(2) - Hnew(2)*constellation_low(int_c_low(2))).^2 + abs(cc(3) - Hnew(3)*constellation_high(int_c_high(3))).^2 + abs(cc(4) - Hnew(4)*constellation_high(int_c_high(4))).^2;
            
            ERRRR(2) = abs(cc(2) - Hnew(2)*constellation_low(int_c_low(2))).^2 + abs(cc(3) - Hnew(3)*constellation_low(int_c_low(3))).^2 + abs(cc(1) - Hnew(1)*constellation_high(int_c_high(1))).^2 +abs(cc(4) - Hnew(4)*constellation_high(int_c_high(4))).^2;
            
            ERRRR(3) = abs(cc(3) - Hnew(3)*constellation_low(int_c_low(3))).^2 + abs(cc(4) - Hnew(4)*constellation_low(int_c_low(4))).^2 + abs(cc(1) - Hnew(1)*constellation_high(int_c_high(1))).^2 +abs(cc(2) - Hnew(2)*constellation_high(int_c_high(2))).^2;
            
            ERRRR(4) = abs(cc(1) - Hnew(1)*constellation_low(int_c_low(1))).^2 + abs(cc(4) - Hnew(4)*constellation_low(int_c_low(4))).^2 + abs(cc(2) - Hnew(2)*constellation_high(int_c_high(2))).^2 + abs(cc(3) - Hnew(3)*constellation_high(int_c_high(3))).^2;
            
            [b comb(k)] = min(ERRRR);
        end
        clear k index_rec sub_car_low_rec sub_car_high_rec  c_low_1 c_high_1 connum_low_rec connum_high_rec;
        for k=1:G
            index_rec = de2bi(comb(k)-1,2);
            sub_car_low_rec=LT(comb(k),1:2);
            sub_car_high_rec=LT(comb(k),3:4);
            c_low_1=c_low((k-1)*nn+1:nn*k);
            c_high_1=c_high((k-1)*nn+1:nn*k);
            connum_low_rec=[nt(c_low_1( sub_car_low_rec(1)));nt(c_low_1( sub_car_low_rec(2)))];
            connum_high_rec=[lt(c_high_1( sub_car_high_rec(1)));lt(c_high_1( sub_car_high_rec(2)))];
            vector_add =[vector_add; connum_low_rec; connum_high_rec];
            demod1 = zeros(2,4);
            demod1(:,1) = step(modu3,connum_low_rec(1,1))';
            demod1(:,2) = step(modu3,connum_low_rec(2,1))';
            demod1(:,3) = step(modu4,connum_high_rec(1,1))';
            demod1(:,4) = step(modu4,connum_high_rec(2,1))';
            
            demod(:,k)=[index_rec.';demod1(:)];
            
        end
        demodf(:,ii)=reshape(demod,p*G,1);
    end
    
    demodf1=reshape(demodf,p*no_symbols*G,1);
    r(c) = nnz(input-demodf1)/(p*no_symbols*G);
end
figure
semilogy([0:1:25],r,'-ok');
grid on;
title('DM-OFDM-MCC');
ylabel('BER');
xlabel('SNR [dB]');
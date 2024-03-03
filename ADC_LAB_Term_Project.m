clc;
clear;
close all;
for z=1:1:5
dict={'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';
    'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z';'a';'b';'c';'d';'e';
    'f';'g';'h';'i';'j';'k';'l';'m';'n';'o';'p';'q';'r';'s';'t';'u';'v';
    'w';'x';'y';'z'};
input=load('sample3text');
x=char(input.test);          %%%%Import the .mat file containing characters%%%%%



snr=-5:1:15;
for o=1:1:length(snr)
clear codes code_bits cw final_txed_bits txed_seq pdesc p i j k l
% % % calculate probability of occurrences of the different characters%%%
char_len=length(x);
p=zeros(1,length(dict));
countchar=0;
for j=1:1:length(dict)
    countchar=0;
    for i=1:1:char_len
    if x(i)==cell2mat(dict(j,1))
     countchar=countchar+1;
    end
    proboccur(j)=countchar/char_len;
    end
end
pdesc=sort(proboccur,'descend');
i=1;
Hs=0;
while (pdesc(i)>0)
    Hs=Hs+pdesc(i)*log2(1/pdesc(i));
    i=i+1;
end
%%
%%-----------------------------SourceCoding---------------------------------%%%%
%Perform fixed length coding of alphabets
codes=randperm(52);
code_bits={};
for i=1:1:length(dict)
code_bits{i,1}=dec2bin(codes(i),6);
end
cw={};
dict=[dict,code_bits];
for i=1:1:char_len
    for j=1:1:length(dict)
        if(x(i)==cell2mat(dict(j,1)))         %Fixed Length Coding length 6 bits
            cw{i,1}=dict(j,2);                
        end
    end
end
final_txed_bits=zeros(char_len,6);
for i=1:1:char_len
    temp=cell2mat(cw{i,1});
        for j=1:1:6 
            if(temp(1,j)=='0')
                final_txed_bits(i,j)=0;      %Fixed Length Bits Generation
            else
                final_txed_bits(i,j)=1;
            end
        end
end

%----------------------Compression(RLC+Huffman)----------------------------------------
txed_seq=final_txed_bits(1,:);
for i=1:1:char_len-1
txed_seq=[txed_seq,final_txed_bits(i+1,:)];      %Final Seq to be transmitted sent to compression
end

%%%%-------------------RLC------------------------%%%%%%%
clear comp_seq count txseq code1
comp_seq=[];
count=[];
j=1;
l=1;
while(j<(length(txed_seq)))                    %Calculate the frequency of 1's and 0's 
i=1;                                           % Take the count and do Huffman coding for the respective count numbers
k=1;
count0=0;
count1=0;
if(txed_seq(j)==0)
    i=j;
    while (txed_seq(i)==0)&&(i<(length(txed_seq)))
        count0=count0+1;
        i=i+1;
    end
    comp_seq(l)=0;
    count(l)=count0;
    l=l+1;
    j=i;
else
    k=j;
    while (txed_seq(k)==1)&&(k<(length(txed_seq)))
        count1=count1+1;
        k=k+1;
    end
    comp_seq(l)=1;                    %comp_seq variable is stored and retrieved at the receiver side 
    count(l)=count1;
    l=l+1;
    j=k;
end
end
count(1,length(count))=count(1,length(count))+1; %Count holds values of decimal numbers depicting frequencies of repetition
LUT=[comp_seq.',count.'];                   %Perform Huffman Coding of the received counts of 1's and 0's
[count_sort,indices]=sort(count,'descend');
temp=count_sort(1);
freq=0;
j=1;
k=1;
l=1;
cb=[];
while(j<length(count_sort))
    k=j;
    while(count_sort(k)==temp)&&(k<length(count_sort))
        freq=freq+1;
        k=k+1;
    end
    cb(l,1)=count_sort(j);
    cb(l,2)=freq/length(count_sort);
    freq=0;
    l=l+1;
    j=k;
    temp=count_sort(k);
end
%%%%%%%%%%%%%%-----Huffmann coding-------------------------%%%%%%%%%%%%%%
clear s b pm ind code1 txseq code p trace i1 i2 I I1 I2 I3 ind prob indice A1
[prob,indice]=sort(cb(:,2).','descend');
s=cb(indice,1).';
A1=count;
 i=1;
p=sort(prob,'descend');
code=[];
b=length(p);
lp=[b];
pm(i,:)=p;
ind(i,:)=1:1:length(s);
code1={};
txseq=[];

%BUILDING PROBABILITY MATRIX
while(length(p)>2) 
 sum1=p(length(p))+p(length(p)-1);
 p=[p(1:length(p)-2),sum1];
 [p,I]=sort(p,'descend');
 i=i+1;
 pm(i,:)=[p,zeros(1,b-length(p))]; %PROBABILITIES MATRIX
 ind(i,:)=[I,zeros(1,b-length(p))];%ARRAY CONTAINING INDICES OF THE SORTED PROBABILITIES
 lp=[lp,length(p)]; %CONTAINS NUMBER OF PROBABILITIES LEFT AT EVERY STEP
end

%TRACING THE PROBABILITIES AND CONSTRUCTING CODEBOOK
c=1;
while(c<=length(s))
    i1=c;
    for i=1:length(s)-1
        trace=ind(i,:);
        i2=find(trace==i1);          %TRACING THE PATH OF PROBABILITY
        I2(i)=i2;
        i1=i2;
        if(i2==lp(i))
            i1=lp(i)-1;
        end
    end
   I3(c,:)=I2;
   c=c+1;
end
for ifinal=1:length(s)
     check1=I3(ifinal,:);
     code=[];
    for i=1:length(s)-1
    
        if(check1(i)==lp(i))
           code=[0 code];                   %ASSIGNING '1's AND '0's
        end
        if(check1(i)==lp(i)-1)
           code=[1 code];
        end
    end
    code1(ifinal,:) = {code};
end
%ENCODING THE FILE USING CODEBOOK
t=A1;
for i=1:length(t)
    for j=1:length(s)
        if(t(i)==s(j))
            temptx=cell2mat(code1(j,:));
            txseq=[txseq temptx] ;
        end
    end
end
compressionratio(z,o)=size(txed_seq,2)/size(txseq,2);            %Compression ratio is calculated
compressedspace_percentage(z,o)=(1-(size(txseq,2)/size(txed_seq,2)))*100;    %Compressed space percentage is calculated
%%        
%--------------------------Channel Encoder---------------------------------
%--------------------------------------------------------------------------

clear n k N D D1 D2 C fbits_txed I P
n=10;
k=5;
D=txseq(1,1:k);
N=floor(length(txseq)/k);
for i=1:1:N-1
D=[D;txseq(1,(i*k)+1:(i*k)+k)];      %Final Seq to be transmitted sent to compression
end
rem=mod(length(txseq),k);
D1=[zeros(1,k-rem),txseq(1,length(txseq)-rem+1:length(txseq))];
D2=[D;D1];
  P=[1,0,1,0,0;0,1,1,0,0;1,1,0,0,0;1,0,0,1,0;0,1,0,1,1];%parity matrx
  I=[1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
  G=[I,P];%Generator matrix
  H=[P.',I];
  C=D2*G;%codevectors
  C_f=mod(C,2);
fbits_txed=C_f(1,:);
for i=1:1:length(C_f)-1
fbits_txed=[fbits_txed,C_f(i+1,:)];      %Final Seq to be transmitted sent to compression
end
%%
%---------------------LC-Mod-AWGN-Demod-Detector----------------------------------------------------
clear m_est N prob1 s A s1 s2 Tb T t sn
%polar NRZ or Antipodal signalling
N=length(fbits_txed);
prob1=nnz(fbits_txed)/N;
A=0.5;
Tb=1;
T=0:0.1:(N*Tb)-0.1;
t=0:0.1:Tb-0.1;
for i=1:1:size(t,2)
s1(i)=A;
end
for i=1:1:size(t,2)
s2(i)=-A;
end
Eb=0;
for i=1:1:size(t,2)
Eb=Eb+s1(i).^2;
end
aphi=[];
aphi=s1./sqrt(Eb/size(t,2));
s=zeros(1,size(t,2));
for k=1:1:length(fbits_txed)
    if (fbits_txed(k)==1)
    s=[s,s1];
    else
    s=[s,s2];
    end
end
figure(1)
plot(T,s(1,size(t,2)+1:size(s,2)));
xlabel('sec');
ylabel('Amplitude');
title('Rectangular impulses based on the generated bitstream');
ylim([-2,2]);
% Adding WGN
sn=awgn(s,snr(o),'measured');
No=[];
No(o)=2*(Eb/size(t,2))*(10^(-snr(o)/10));
figure(2)
plot(T,sn(1,size(t,2)+1:size(s,2)));
xlabel('sec');
ylabel('Final signal amplitude with noise');
title('signal + WGN');
ylim([-10,10]);
% Correlator bank
%%%%%Without AWGN noise%%%%%%%
% y=[];
% y=s(1,size(t,2)+1:size(s,2));
% s_rec=[];
% for k=1:1:N
% s_rec=y(1,(size(t,2)*(k-1))+1:size(t,2)*k);
% x=sum(s_rec.*aphi);
% if (round(x-sqrt(Eb))==0)
%     m_est(k)=1;
% else
%     m_est(k)=0;
% end
% end
% pe=nnz(bits_gen-m_est)/N;
%%%%%With AWGN Noise%%%%%%%
y=sn(1,size(t,2)+1:size(s,2));
alpha(o)=(No(o)/(4*sqrt(Eb/size(t,2))))*log((1-prob1)/prob1);
for k=1:1:N
    s_rec=[];
    s_rec=y(1,(size(t,2)*(k-1))+1:size(t,2)*k);
   y_rec=conv(s_rec,aphi);
    r1=y_rec(size(t,2)+1)./size(t,2);
    if (r1>alpha(o))
        m_est(k)=1;
    elseif(r1<alpha(o))
        m_est(k)=0;
    end
end
%%
%-----------------------------Channel Decoder-------------------------------------------------------
%---------------------------------------------------------------------------------------------------
clear R N2 bits_corr fbit_rxed S totalbit_rxed totalbit_rxedf
n=10;
k=5;
R=m_est(1,1:n);
N2=length(m_est)/n;           %block size
for i=1:1:N2-1
R=[R;m_est(1,(i*n)+1:(i*n)+n)];      %received vector R
end
 bits_corr=zeros(N2,n);
  fbit_rxed=zeros(N2,k);
  S=mod(R*(H.'),2);                  %finding the syndrome
  Ht=H.';
  j=1;
  for i=1:1:N2
      if (S(i,1)+S(i,2)+S(i,3)+S(i,4)+S(i,5))==0
          bits_corr(i,:)=R(i,:);
      else
          ht=size(H.');
          for l=1:1:ht(1)
              if S(i,:)==Ht(l,:)
                  R(i,l)=not(R(i,l));
                  break;
              else
                  R(i,l)=R(i,l);
                  continue;
              end
          end
          bits_corr(i,:)=R(i,:);
      end
  end
  fbit_rxed=bits_corr(:,1:k);
  totalbit_rxed=fbit_rxed(1,:);
  for i=1:1:N2-1
totalbit_rxed=[totalbit_rxed,fbit_rxed(i+1,:)];  
  end
  totalbit_rxedf=[totalbit_rxed(1,1:length(totalbit_rxed)-k),totalbit_rxed(1,length(totalbit_rxed)-rem+1:length(totalbit_rxed))];
%%
%----------------------------------Source Decoder----------------------------------------------------  
%----------------------------------------------------------------------------------------------------
%%%%%%%%Huffman Decoder%%%%%%%%%%%%%%%%
clear dec decoded_seq rhuf
it=1;
dec=[];
s1=cb(indice,1).';
min5=length(code1{1,1});
rhuf=totalbit_rxedf;
for i=1:length(s1)
        comp1=code1{i,1};
        if(length(comp1)<min5)
            min5=length(comp1);
        end
end
while(it<=length(rhuf))
    for i=1:length(s1)
       temptx=cell2mat(code1(i,:)); %TAKING CODEWORD OF EACH SYMBOl
        if((it+length(temptx)-1)<=length(rhuf))
           comp=rhuf(it:it+length(temptx)-1); %TAKING BITS FROM RECEIVED SEQUENCE EQUAL TO NUMBER OF BITS OF THE CODEWORD THAT WE ARE COMPARING TO
          comp2(i)=sum(xor(comp,temptx));  %CALCULATING HAMMING DISTANCE
       end
    end
    min1=min(comp2);                      %TAKING MINIMUM HAMMING DISTANCE
    mm=find(comp2==min1);                 %DECODING THE SYMBOL WICH GAVE MINIMUM HAMMING DISTANCE
    it=it+length(cell2mat(code1(mm(1),:)));
    dec=[dec s1(mm)];                       %DECODED FILE
       if (length(rhuf)-it<=min5-1)
           it=length(rhuf)+1;
           break;
       end
end
%%%%%%%%RlC Decoding%%%%%%%%%%%%%%%%%
decoded_seq=[];
for m=1:1:length(comp_seq)
    t=1;
     while (t<=dec(m))
        decoded_seq=[decoded_seq,comp_seq(m)]; %%%%%stored comp_seq is used here
        t=t+1;
    end
    if(m==length(dec))
        break;
    else
        continue;
    end
    
   
end

if(length(decoded_seq)<length(txed_seq))
    seq_final=[decoded_seq,ones(1,(size(txed_seq,2)-size(decoded_seq,2)))];
elseif(length(decoded_seq)>length(txed_seq))
    seq_final=decoded_seq(1,1:length(txed_seq));
else
    seq_final=decoded_seq;
end
l=1;
%%
%%%%%%%%%Fixed length decoding%%%%%%%%%%%%%%%%%%%%%%%%%%
clear recchar
for j=1:6:length(seq_final)
    k=1;
    tempdec=(2^5)*seq_final(1,j)+(2^4)*seq_final(1,j+1)+(2^3)*seq_final(1,j+2)+(2^2)*seq_final(1,j+3)+(2^1)*seq_final(1,j+4)+seq_final(1,j+5);
    while (k<=length(dict))
        if (bin2dec(cell2mat(dict(k,2)))==tempdec)
            recchar(1,l)=char(dict(k,1));    %%%%%Retrieved Characters %%%%%%%%%
            l=l+1;
            k=k+1;
        else
            k=k+1;
        end
    end
end
%%
%%%%%%%%%%%%%Symbol Error rate calculation%%%%%%%%%%%%%%%%%%%
no_errors(z,o)=0;
for q=1:1:length(recchar)
    if(x(q)==recchar(q))
        no_errors(z,o)=no_errors(z,o);
    else
        no_errors(z,o)=no_errors(z,o)+1;
    end
end
        
end

end
for w=1:1:length(snr)
SER(w)=sum(no_errors(:,w).')./(char_len*5);
end
%%
figure
semilogy(snr,SER)
xlabel('snr')
ylabel('SER');
title('Symbol error rate');
grid on;
disp(compressionratio);
  
rng('default')
clear
clc
%% channel generation
nA  = 4;
nB = 4;
nE = 4;

maxIter = 100;

SNRdB = 15;
P0 = 10^(SNRdB/10);
PAPC = ones(nA,1)*P0/nA*1.2;
Titer = 100;
toatltime=0;
toatltimePABR=0;

APGmat = zeros(Titer,1);
Matdat=[];
r = 0; % transmit antenna correlation
phib = 0;
phie = 90/180*(pi);
Rb = toeplitz((r*exp(1i*phib)).^(0:nA-1));
Re = toeplitz((r*exp(1i*phie)).^(0:nA-1));

mygamma = 0.1; % this number is to reduce the strength of He channel to obtain more degraded channels


nChannels = 100; % total number of generated channels

errorseq = zeros(maxIter,nChannels);
iNonDegradChan = 0;
for iChan =1:nChannels
    Hb =(randn(nB,nA)+1i*randn(nB,nA))/sqrt(2)*sqrtm(Rb);
    He = sqrt(mygamma)*(randn(nE,nA)+1i*randn(nE,nA))/sqrt(2)*sqrtm(Re);
    H  = [Hb;He];
    Delta = (Hb'*Hb-He'*He);
    if(min(real(eig(Delta)))<0) % if the channel is nondegraded
        iNonDegradChan = iNonDegradChan +1;
        [SecrecyCapAdaptiveMomentum,objseq_SecrecyCapAdaptiveMomentum] = Algorithm1(Hb,He,nA,nB,nE,maxIter,P0);
        [SecrecyCap] = SecrecyCapacityPBRA(H,Hb,He,nA,nB,nE,maxIter,P0);
        errorseq(1:length(objseq_SecrecyCapAdaptiveMomentum),iNonDegradChan)= SecrecyCap-objseq_SecrecyCapAdaptiveMomentum;
        
    end
end
fprintf('Total of generated channels: %d\n',nChannels)
fprintf('Number of nondegraded channels: %d\n',iNonDegradChan)
errorseq(:,iNonDegradChan+1:end)=[];
meanerrorseq = mean(errorseq,2);
semilogy(errorseq,'LineWidth',.5,'Color','#D3D3D3')
hold on
semilogy(meanerrorseq,'r','LineWidth',2)
grid on
ylim([9e-6,10])
xlabel('Iteration count')
ylabel('Optimality gap')
saveas(gcf, '../../results/Fig1b.png')
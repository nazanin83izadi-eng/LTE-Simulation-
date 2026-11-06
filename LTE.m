%% LTE SIMULATION
enb.NDLRB = 25 ;       %5MHZ BW
enb.CellRefP=1;        % 1 antenna port (1 2 4)
enb.NCellID = 0;       % AZ 0 TA 503
enb.CyclicPrefix = 'Normal';   % cp= extended/normal
enb.CFI = 2;  % 2 symbol allocated to control channel
enb.Ng = 'sixth';
enb.DuplexMode = 'FDD';    
%%enb.NSubframe = 0;      % subframe number = 0,....,9

Imcs = 15;
[pdsch.Modulation,Qm , itbs] =  mcsModulation(Imcs);

%%PDSCH PARAMETERS
pdsch.RNTI = 1;
pdsch.PRBSet = (0 : enb.NDLRB-1).' ;   % PHY RESOURSE BLOCK        ./:ستونی
pdsch.TxScheme = 'Port0' ;
pdsch.NLayers = 1;

%%PDCCH PARAMETERS
pdcch.NDLRB = enb.NDLRB;
pdcch.RNTI = 1;
pdcch.PDCCHFormat =0;


%%SNR (AWGN)
SNRdBrange = 0:1:20;      % اثر روی نویز کانال 
BER = zeros(size(SNRdBrange));
numFrames = 10;

 for SNR = 1:length(SNRdBrange)

     SNRdB = SNRdBrange(SNR);

     totalErrBits = 0;
     totalBits    = 0;

for frame = 1:numFrames

fullGrid = [];
origBits = cell(10,1);             % for compare subframe in Rx
RxBits   = cell(10,1);       
crcErr   = zeros(1,10);

for sf = 0:9
    enb.NSubframe  = sf;

    INDpdsch = ltePDSCHIndices(enb , pdsch ,pdsch.PRBSet );

nRE= size(INDpdsch , 1);     % number of resource elements
bitlength = nRE * Qm;

%%PDCCH
INDpdcch = ltePDCCHIndices(enb);
npdcch = length(INDpdcch); % number of elements
npdcchbit = npdcch *2 ; % QPSK
pdcchBit = randi ([0 1] ,npdcchbit ,1);

origPDCCHBits{sf+1} = pdcchBit;

%%TransportBLK
nprb = enb.NDLRB;                    % 1 user
tbs = lteTBS(nprb , itbs);

%%CRC
blk = randi ([0 1] , tbs ,1);  %making a block
origBits{sf+1} = blk;
poly = '24A';
blkcrc = lteCRCEncode(blk , poly);

%%segmentation
cbs = lteCodeBlockSegment(blkcrc);

%%Turbo code
turbocode = lteTurboEncode(cbs);

%%Matching Rate
rv = 0;                  %% redundancy version
ratematch = lteRateMatchTurbo(turbocode , bitlength, rv);


%%SYMBOL
symbols = ltePDSCH(enb , pdsch , ratematch);
symbolsc = ltePDCCH(enb ,pdcchBit);


%%GRID
Grid = lteResourceGrid(enb);    %2*2 MATRIX (خالی)
Grid(INDpdsch)= symbols;
Grid(INDpdcch)= symbolsc;

%%CRS PILOT
refsignal = lteCellRS(enb);
refsignalInd = lteCellRSIndices(enb);
Grid(refsignalInd) = refsignal;

%%PSS/SSS
if sf == 0 || sf==5
    pss = ltePSS(enb);
    pssInd = ltePSSIndices(enb);
    Grid(pssInd) = pss;

   
    sss= lteSSS(enb);
    sssInd = lteSSSIndices(enb);
    Grid(sssInd) = sss;

end

fullGrid = [fullGrid , Grid];

end



%%OFDM MODULATION
enb.NSubframe = 0;
[waveform,info] = lteOFDMModulate(enb,fullGrid );

%%LTEFADINGCHANNEL
cfg.Seed = frame;
cfg.NRxAnts = 1;
cfg.DelayProfile = 'EPA';     %%EPA/EVA/ETU
cfg.DopplerFreq = 5;      %% based on users speed
cfg.InitTime = 0;
cfg.MIMOCorrelation = 'Low'; 
cfg.NTerms = 16;
cfg.ModelType = 'GMEDS'; 
cfg.InitPhase = 'Random';
cfg.NormalizePathGains = 'On';
cfg.NormalizeTxAnts = 'On';

cfg.SamplingRate = info.SamplingRate;  
%%disp(cfg);
fading = lteFadingChannel(cfg,waveform);

%AWGN CHANNEL
N0 = 1/(sqrt(2.0*double(enb.CellRefP)*double(info.Nfft)) * (10^(SNRdB/20)));
    noise = N0 * complex(randn(size(fading)), randn(size(fading)));


    Rxwaveform = fading + noise;



 %RECIEVER
 RxGrid = lteOFDMDemodulate(enb ,Rxwaveform );


 %channel estimate   *****
cec.PilotAverage = 'UserDefined'; % Pilot averaging method
cec.FreqWindow = 19;               % Frequency averaging window in REs
cec.TimeWindow = 5;               % Time averaging window in REs
cec.InterpType = 'linear';         % linear interpolation
cec.InterpWinSize = 1;            
cec.InterpWindow = 'Centred';     % Interpolation windowing method

%%decoding

Rxbit = cell(10 ,1);
RxPDCCHbit = cell(10,1);


for sf = 0:9
    enb.NSubframe = sf;


    colS = sf*14 + 1;
    colE = sf*14 + 14;
    sfRxGrid = RxGrid(: , colS:colE , :);

  [Channelest, noiseEst] = lteDLChannelEstimate(enb,cec,sfRxGrid);  


%%pdcch   ****
pdcchIND =ltePDCCHIndices(enb);

[pdcchRx , pdcchest] = lteExtractResources(pdcchIND , sfRxGrid , Channelest);
pdccheq = lteEqualizeMMSE(pdcchRx ,pdcchest,noiseEst );
%چک بشهه ****************
pdcchSOFTbit = ltePDCCHDecode(enb , pdccheq);

pdcchbits{sf + 1} = pdcchSOFTbit;



%%PDSCH

INDpdsch = ltePDSCHIndices(enb,pdsch,pdsch.PRBSet);

nRE = size(INDpdsch,1);
bitlength = nRE*Qm;

[pdschRx,pdschest] = lteExtractResources(INDpdsch,sfRxGrid,Channelest);

pdscheq = lteEqualizeMMSE(pdschRx,pdschest,noiseEst); 


pdschSoftBits = ltePDSCHDecode(enb,pdsch,pdscheq);
pdschSoftBits = pdschSoftBits{1};

rv = 0;

raterec = lteRateRecoverTurbo( pdschSoftBits, tbs, rv);

turboDecode = lteTurboDecode(raterec);

deseg = lteCodeBlockDesegment(turboDecode,tbs+24);

[RxBits{sf+1},crcErr(sf+1)] = lteCRCDecode(deseg,poly);


%%BER
 errBits = sum(RxBits{sf+1} ~= origBits{sf+1});
        totalErrBits = totalErrBits + errBits;
        totalBits    = totalBits + length(origBits{sf+1});
    end

end   

    BER(SNR) = totalErrBits / totalBits;

end



InfS = info.SamplingRate;   

figure;
pspectrum(waveform, InfS);
title('(Tx Waveform)');

figure;
pspectrum(fading, InfS); title('After Fading Channel');

figure;
pspectrum(Rxwaveform, InfS);
title('(Rx Waveform + AWGN)');

%%PLOT BER
figure;
semilogy(SNRdBrange, BER, '-o', 'LineWidth', 2);
grid on;
xlabel('SNR');
ylabel('BER');

fprintf( 'Imcs :%d\n',Imcs);
fprintf( 'Qm :%d\n',Qm);
disp(cfg);
disp(cec);

%%Function of MCS
function [modulation, Qm , itbs] = mcsModulation(Imcs)
 if Imcs >=0 && Imcs<= 9
     modulation = 'QPSK';
     Qm = 2;
     itbs = Imcs;
 elseif Imcs >=10 && Imcs <= 16
     modulation = '16QAM';
     Qm = 4;
     itbs = Imcs-1;
 elseif Imcs >=17 && Imcs<= 28
     modulation = '64QAM';
     Qm = 6;
     itbs = Imcs-2;
 end
end





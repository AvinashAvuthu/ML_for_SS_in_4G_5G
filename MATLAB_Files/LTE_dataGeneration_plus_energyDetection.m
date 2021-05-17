clear;
clc;
addpath('functions')

ifWrite2File = 1;                                                           % 1 if data is to be written to file, 0 otherwise
numfiles = 2;                                                               % 2 files - first for ML trainig, second for ML testing

numUsers = 20;                                                            % number of users/number of different applications for allocation

polar = 0;                                                                  % bits polarity
modType = 'QPSK';                                                           % type of modulation used
numBitsPerSymb = NumberBits(modType);                                       % number of bits per modulation symbol
                                                      

% OFDM variables
FFTLength = 1024;                                                           % numberc of fft samples used in ofdm
numSubcarriers = 12;                                                        % number of subcarriers per LTE resource block
numCPSamples = 72;                                                       % number of samples of OFDM symbol used for cyclic prefix
timeSlotsNum = 80;                                                          % total number of time slots used in program
resBlocksNum = 50;                                                          % number of resource blocks for chosen band (for example resBlocksNum = 50 for 10 MHz LTE channel)
ofdmSymbPerTimeSlot = 1;                                             % number of OFDM sybmols per one time slot

resourceAlloc = zeros(resBlocksNum+1, timeSlotsNum+1);                      % table of resources assigned to users (for visualization purposes)
resourceAllocShort = zeros(numUsers,4);                                     % table of resources assigned to users (for chcecking spectrum occupancy)
                                                                            %   columns: first freq, last freq, first time slot, last time slot
                                                                            %   rows: users

detectedRes = zeros(resBlocksNum*numSubcarriers, timeSlotsNum);    % making a 600*80 matrix
allocPossible = zeros(1, numUsers);                                          % determines if resources allocation is possible at a given time slot

LTETransmData = zeros(1,timeSlotsNum*ofdmSymbPerTimeSlot...     % making a 1*93440 matrix total samples in time domain
    *(FFTLength+numCPSamples));                                              % LTE data transmitted in time domain

Pf = 0.1;                                                                    % Probability of false alarm - used in ED threshold calculation

%% channel parameters

% SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR = -40:1:20;

% Fading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 14.336e6;            %old value: 14.336e6                                           % Sample rate (Hz)
pathDelays = [0 30e-9 70e-9 90e-9 110e-9 190e-9 410e-9];            % Path delays (s)
pathPower = [0 -1 -2 -3 -8 -17.2 -20.8];                            % Path power (dB)
fD = 0;                                                             % Maximum Doppler shift (Hz)
numSamples = length(LTETransmData);                                 % Number of samples


%% signal generation

for numsamples=1:41
    disp("Number of iteration:")
    disp(numsamples)
    for fileIdx=1:numfiles
        LTETransmData = zeros(1,timeSlotsNum*ofdmSymbPerTimeSlot*...
            (FFTLength+numCPSamples));
        resourceAlloc = zeros(resBlocksNum+1, timeSlotsNum+1);
        resourceAllocShort = zeros(numUsers,4);
        allocPossible = zeros(1, numUsers);
        detect2 = zeros(1, length(SNR));
        falseAl2 = zeros(1, length(SNR));

        % fading channel impulse response generation
        rchan = comm.RayleighChannel('SampleRate',fs, ...
            'PathDelays',pathDelays,'AveragePathGains',pathPower, ...
            'MaximumDopplerShift',fD,'FadingTechnique','Sum of sinusoids');

        disp('   resource allocation');

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LTE-RB generation
        for i = 1:numUsers
            while allocPossible(i)== 0
                numTimeSlots = randomParam(4, 30, true, 1, timeSlotsNum);                      % number of OFDM symbols (time domain)
                firstTimeSlot = randomParam(40, 40, true, 1, timeSlotsNum-numTimeSlots+1);     % first time slot index
                numResourceBlocks = randomParam(4, 30, true, 1, resBlocksNum);               % number of resource blocks (frequency domain)
                firstRB = randomParam(25, 20, true, 1, resBlocksNum-numResourceBlocks+1);    % first resource block index
                lastTimeSlot = firstTimeSlot + numTimeSlots - 1;                               % last time slot index
                lastRB = firstRB + numResourceBlocks - 1;                                    % last resource block index
                numBits = ofdmSymbPerTimeSlot*numSubcarriers*numBitsPerSymb*...
                    numTimeSlots*numResourceBlocks;                                            % number of data bits to generate in the beginning of program

                [resourceAllocShort, allocPossible(i)] = CheckResourcesAllocation(firstRB, ...
                    lastRB, firstTimeSlot, lastTimeSlot, resourceAllocShort, ...
                    i, timeSlotsNum, resBlocksNum);
            end

            if allocPossible(i)
                resourceAlloc(firstRB:lastRB, firstTimeSlot:lastTimeSlot) = ...
                    i*ones(numResourceBlocks, numTimeSlots);

                % calculating frequency subcarriers indexes for resource blocks allocated for given user
                if firstRB > lastRB
                    disp([firstRB lastRB])
                end
                [firstSubcarrier, lastSubcarrier] = SubcarrierAllocation([firstRB lastRB], ...
                    FFTLength, numSubcarriers, resBlocksNum);

                %transmiting ofdm data
                ofdmData = OneLTEUserTransm(numBits, polar, modType, ...
                    numTimeSlots, FFTLength, numCPSamples, firstSubcarrier,...
                    ofdmSymbPerTimeSlot);
                ofdmVectL = numTimeSlots*ofdmSymbPerTimeSlot*length(ofdmData(:,1));
                ofdmVect = reshape(ofdmData, [1,ofdmVectL]);
                LTETransmData = transmissionLTE(ofdmVect, firstTimeSlot, numTimeSlots, ...
                    LTETransmData, FFTLength, numCPSamples, ofdmSymbPerTimeSlot);
            end

        end

        %  ALLOCATED RESOURCES VISUALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                figure('Renderer', 'painters', 'Position', [200 200 320 200])
                s = surf(sign(resourceAlloc));
                view(0, 90)
                s.EdgeColor = 'none';
                xlabel('time')
                ylabel('frequency')
                xlim([1,timeSlotsNum+1])
                ylim([1,resBlocksNum+1])
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % even out power
        [LTETransmData, powers] = rescaleLTEData(LTETransmData, resourceAlloc,...
            FFTLength, numCPSamples, ofdmSymbPerTimeSlot);

        %% AWGN channel

        for SNRind = 1:length(SNR)
            disp(['      SNR = ' num2str(SNR(SNRind))])
            variance = 0.002;

            %% transmition through channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % fading
            LTETransmDataChannel = rchan(LTETransmData.').';
            [C, LTETransmDataChannel] = changeSigPow(variance, ...
                SNR(SNRind), LTETransmDataChannel);
            noise =  sqrt(variance/2)*(randn(size(LTETransmDataChannel)) ...
                + 1i*randn(size(LTETransmDataChannel)));
            noisySignal = LTETransmDataChannel + noise;
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Energy Detection

            [detectedRes, energyMtrx] = energyDetection(FFTLength, ...
                numSubcarriers, resBlocksNum, ...
                timeSlotsNum, ofdmSymbPerTimeSlot, numCPSamples, Pf, noisySignal, ...
                variance);
            decisionTable = detectedRes;
            [detectProc, falseAlarmProc] = ...
                detectionCheck(sign(resourceAlloc(1:resBlocksNum, 1:timeSlotsNum)), ...
                detectedRes);

            detect2(SNRind) = detect2(SNRind) + detectProc;
            falseAl2(SNRind) = falseAl2(SNRind) + falseAlarmProc;

            %%  DETECTED RESOURCES VISUALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         figure('Renderer', 'painters', 'Position', [200 200 320 200])
%                         s2 = surf(decisionTable);
%                         s2.EdgeColor = 'none';
%                         view(0, 90)
%                         xlabel('time')
%                         ylabel('frequency')
%                         zlabel('user')
%                         xlim([1,timeSlotsNum+1])
%                         ylim([1,resBlocksNum+1])
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Write to file
            if ifWrite2File == 1
                if mod(fileIdx,2) == 1
                    fileEnd = '.txt';
                else
                    fileEnd = '_test.txt';
                end
                if SNR(SNRind)<0
                    fileName = ['data_' num2str(abs(SNR(SNRind))) fileEnd];
                else
                    fileName = ['data' num2str(abs(SNR(SNRind))) fileEnd];
                end
                fileID = fopen(fileName, 'a+');

                rB = 1;
                tS = 1;
                [maxRB, maxTS] = size(decisionTable);

                forget_Factor = 0.9;

                for rB = 1:resBlocksNum
                    ff_value = 0;

                    for tS = 1:timeSlotsNum

                        try
                            neighbourUp1 = decisionTable(rB+1, tS);
                        catch err
                            if err.identifier == "MATLAB:badsubscript"
                                neighbourUp1 = 0;
                            else
                                disp(err.identifier)
                            end
                        end
                        try
                            neighbourDown1 = decisionTable(rB-1, tS);
                        catch err
                            if err.identifier == "MATLAB:badsubscript"
                                neighbourDown1 = 0;
                            else
                                disp(err.identifier)
                            end
                        end

                        wallNeighbours1 = neighbourUp1 + neighbourDown1;

                        try
                            cornerUp = decisionTable(rB+1, tS-1);
                        catch err
                            if err.identifier == "MATLAB:badsubscript"
                                cornerUp = 0;
                            else
                                disp(err.identifier)
                            end
                        end
                        try
                            cornerDown= decisionTable(rB-1, tS-1);
                        catch err
                            if err.identifier == "MATLAB:badsubscript"
                                cornerDown= 0;
                            else
                                disp(err.identifier)
                            end
                        end

                        cornerNeighbours = cornerUp + cornerDown;


                        ff_value = decisionTable(rB, tS) + forget_Factor*ff_value;



                        fprintf(fileID, '%d %d %d %d %d %d %f\n', ...
                            rB, tS, sign(resourceAlloc(rB,tS)), decisionTable(rB, tS), ...
                            wallNeighbours1, ...
                            cornerNeighbours, ff_value);
                    end
                end

                fclose(fileID);
            end

        end

    % plot enerdy detecion probability of detection and probability of false alarm
        figure
        plot(SNR, detect2, 'b', 'LineWidth', 1)
        hold on
        plot(SNR, falseAl2, 'r', 'LineWidth', 1)
        xlabel('SNR [dB]')
        ylabel('Probability [%]')
        legend('probability of detection', 'probability of false alarm');
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    end
end
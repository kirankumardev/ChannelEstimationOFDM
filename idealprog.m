close all 
clear all

% Signal-to-noise ratios (SNRs) to be simulated
% SNR in dB
SNRindB = [0:1:20]; 

% SNR linear
SNR = 10.^(SNRindB/10);  

% Variance of additive white Gaussian noise (AWGN) at receiver
sigma2_n = 1./SNR;              


% Channel model parameters

% Number of taps within the channel impulse response
K = 2; 
% Factor for exponentially decaying power profile: 
c_att = 1; 
 % Number of channel realizations to be simulated
Nreal = 10^6;  


% OFDM parameters
% Number of carriers/ FFT-size 
Nc = 128; 

% Code rate of employed repetition code                    
R = 1/4;

% intlv = 1 for evenly spaced pilots
% intlv = 0 for closely spaced pilots
intlv = 1;          
                    
% Generate index vector for interleaving with maximum distance pattern for
% interleaved pilots
if ( intlv )&( R < 1 )
    index_matrix = [];
    for kk=1:1/R
        index_matrix = [index_matrix [(kk-1)*Nc*R+1:kk*Nc*R]' ];
    end
    index = reshape(index_matrix',1,Nc);
elseif ( ~intlv )&( R < 1 )
    index_matrix = reshape([1:Nc],1/R,Nc*R)';
end

% Calculate variances of channel taps according to exponetial decaying power profile
var_ch = exp(-[0:K-1]); 

% Normalize overall average channel power to one
var_ch = var_ch/sum(var_ch);       

% Variable for calculating Bit Error Rate
err_count = zeros(1,length(SNRindB));

% Simulation for BER
for ii = 1:length(SNRindB)
    ii
    for jj = 1:Nreal
 
        %Transmitter 
        % Generate random data vector of length Nc*R with elements {-1,+1}
        U = 2*round(rand(1,R*Nc))-1;
         
        % Perform repetition encoding [+1 -1 ...] --> [+1 +1 -1 -1 ...]
        if R < 1
            X = kron(U,ones(1,1/R));
        else
            X = U;
        end
 
        % Perform interleaving   
        if ( intlv )&( R < 1 )
            X(index) = X;
        end             

       % IFFT of current OFDM symbol X including normalization factor(Unitary IFFT) 
        x = ifft(X)*sqrt(Nc);

        % Add cyclic prefix of length Nch-1(Transmitted sequence)
        x = [ x(end-K+2:end) x ];

       
        % Channel part
        % Generate random channel coefficients (Rayleigh fading)
        h = sqrt(0.5)*( randn(1,K) + j*randn(1,K) ) .* sqrt(var_ch); 

        % Calculate corresponding frequency response needed for receiver part
        % Zero-padded channel impulse response (length Nc)
        h_zp = [h zeros(1,Nc-K)]; 
        
        % Corresponding FFT 
        H = fft(h_zp);                               

        % Received sequence obtained through Convolution with channel impulse response
        y = conv(x,h);

        % Add AWGN 
        n = sqrt(0.5*sigma2_n(ii)) * ( randn(1,length(y)) + j*randn(1,length(y)) );
        y = y + n;

        % Discard last Nch-1 received values resulting from convolution
        y(end-K+2:end) = [];   

        % Receiver
        % Remove cyclic prefix
        y(1:K-1) = [];   
 
        % FFT of received vector including normalization factor (--> unitary FFT)
        Y = fft(y)/sqrt(Nc);

        % Perform deinterleaving and repetition decoding (--> maximum ratio combining)
        % Derotation and weighting
        Z = conj(H) .* Y;                   
        if R == 1
            Uhat = sign(real(Z));
        else
            % Perform deinterleaving
            matrix_help = Z(index_matrix);
            % Maximum ratio combining
            Z_mrc = sum(matrix_help,2); 
            % Hard decision on information symbols
            Uhat = sign(real(Z_mrc))';      
        end 

        % Bit error count
        err_count(ii) =  err_count(ii) + sum(abs(Uhat-U))/2;

    end % loop over all channel realizations
end % loop over all SNR values 

% Calculate final bit error rate (BER)
ber = err_count/(R*Nc*Nreal);


% Plot for knowing the Channel Frequency Response
% Plot one random realization of channel frequency response for last
% simulation
figure
set(gca,'FontSize',12);
h=plot(abs(H),'b-');
set(h,'LineWidth',1);
hold on
if R < 1
    for ii=1:1/R
        if intlv
           h=stem(index(ii),abs(H(index(ii))),'r');
           set(h,'MarkerSize',7,'LineWidth',1);
        else
           h=stem(ii,abs(H(ii)),'r');
           set(h,'MarkerSize',7,'LineWidth',1);
        end
    end
end
xlabel('Carrier No.')
ylabel('Magnitude of Frequency Response of the Channel(H)')
axis([1 Nc 0 1.1*max(abs(H))])
grid on

% BER plot
figure
set(gca,'FontSize',12);
h=semilogy(SNRindB,ber,'ko--');
set(h,'MarkerSize',7,'LineWidth',1);
hold on
legend('Simulation for M = 4 and K = 2')
xlabel('SNR = 1/\sigma_n^2 in dB')
ylabel('BER')
axis([SNRindB(1) SNRindB(end) 0.0001 1])
grid on

 
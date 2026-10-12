% Reference signal the intermodulation script expects but TestScripts.zip does not ship:
% the ideal 44.1 kHz rendering of intermodulation_sine (60 Hz at -6 dBFS + 7 kHz at -18.0412 dBFS).
% The frequencies are the exact values used by GEN_intermodulation_60Hz_and_7kHz.m.
t = 0:1/44100:10;
t = t(1:end-1);
y = 10^(-6/20) * sin(2*pi*64.599609375*t) + 10^(-18.0412/20) * sin(2*pi*6998.291015625*t);
audiowrite('intermodulation_sine-InternalUseAt44kHz.wav', y(:), 44100, 'BitsPerSample', 64);

function H = constellations( D, fs , DRAW)
  % returns an 2-D array with its rows as [f1, f2, dt, t1].
  %   D:  1-D waveform data
  %   fs: sampling rate for D
  %   DRAW: =1: draw spectrum. might be slow!
  
  if nargin == 2; DRAW = 0; end;
  
  %% CONFIGURABLE PARAMETERS %%
  
  % configurable short-time Fourier transform parameters
  STFT_WINDOW  = 256;  % STFT time window. don't know what it is :P
  STFT_OVERLAP = 220;  % # of STFT samples overlapped with last STFT.
  STFT_NSAMPLE = 256;  % # of STFT samples
  
  % peak-marking parameters
  PEAK_NUM = 5; % find largest PEAKNUM peaks each time bin.

  % constellation parameters
  %
  % window   ___________________
  % offset   |        x         |
  %      \   |     x         x  |
  %     x----|         x        |  window height (freq bins)
  %          |   x          x   |
  %          |__________________|
  %         window width (time bins)
  %

  WINDOW_HEIGHT = 10; % unit: frequency bins
  WINDOW_WIDTH = 10;  % unit: time bins
  WINDOW_OFFSET = 1;  % unit: time bins
  
  %% spectrum analysis %%
  
  [S,F,T,P] = spectrogram(D, STFT_WINDOW, STFT_OVERLAP, STFT_NSAMPLE, fs); 
  % S: ignored. 
  % F: freq bin -> real frequency in Hz. Can be used as frequency axis
  % T: time bin -> real time in seconds. Can be used as time axis.
  % P: (freq bin size)x(time bin size) matrix. Its value is power.
  
  fprintf('%d time bins, %d frequency bins.\n', length(T), length(F));
  time_res = T(2)-T(1)
  fprintf('Time resolution = %f sec\n', time_res);
  fprintf('Frequency resolution = %f Hz\n', F(2)-F(1));
  
  P = 10*log10(abs(P)); % work in decibels.

  if DRAW % draw spectrum. might be slow!
    spectrogram(D, STFT_WINDOW, STFT_OVERLAP, STFT_NSAMPLE, fs, 'yaxis'); 
    hold on;
  end
  
  %% find local maxes %%
  
  % local maxes array -- peaks
  % peaks(i, j, 1) = jth peak value of time bin i
  % peaks(i, j, 2) = location (freq bin) of jth peak of time bin i
  peaks = zeros(length(T), PEAK_NUM);
  peaks(:,:,2) = zeros(length(T), PEAK_NUM);
  fprintf('Processing');
  for i = 1:length(T) % for each time bin
    %find local maximals that is 3dB higher than surroundings.
    [pks, loc] = findpeaks(P(:, i),'THRESHOLD' , 3, 'SORTSTR', 'descend');
    
    % only keeps largest PEAK_NUM peaks
    nPeaks = length(pks);
    if nPeaks > PEAK_NUM
      pks = pks(1:PEAK_NUM); loc = loc(1:PEAK_NUM);
      nPeaks = PEAK_NUM;
    end
    peaks(i,1:nPeaks, 1) = pks;
    peaks(i,1:nPeaks, 2) = loc;
      
    % plot peaks
    if DRAW;  stem3(ones(1,length(pks)) .* T(i), F(loc), pks); end;
    
    % progress dots: output a '.' everytime when 100 time bins is processed
    if mod(i, 100) == 0; fprintf('.'); end;
  end
  fprintf('\n');

  %disp(peaks);
  
  %% packing local maxes into H %%
  
  fprintf('Generating Pairs');
  peaks = peaks(:,:,2); % discard peak value, only loc remains
  
  prealloc_size = PEAK_NUM * length(T);
  H = zeros(prealloc_size , 4);
  % each row are [f1, f2, dt, t1]
  
  Hptr = 1; % current pointer to the end of H
  
  for t = 1:length(T)-(WINDOW_OFFSET + WINDOW_WIDTH)
    t_min = t + WINDOW_OFFSET;
    t_max = t + WINDOW_OFFSET + WINDOW_WIDTH;    
    pks = peaks(t, peaks(t, :) > 0 ); % peaks in this time bin (t)
    target_peaks = peaks(t_min:t_max, :); % peaks in time bin t_min~t_max
    
    for f1 = pks % for each peaks in this time bin
      f_min = f1 - WINDOW_HEIGHT / 2;
      f_max = f1 + WINDOW_HEIGHT / 2;
      
      f2_index = target_peaks > f_min & target_peaks < f_max; % logical arr
      f2 = target_peaks(f2_index);
      dt = mod(find(f2_index)-1, size(target_peaks,1)) + WINDOW_OFFSET;  % row num of logical arr = dt arr
      n = length(f2); % size of f2 = size of dt
      
      H(Hptr:Hptr+n-1, :) = [ones(n,1).*f1,f2,dt,ones(n,1).*t];
      Hptr = Hptr + n;
    end
    
    if mod(t, 100) == 0; fprintf('.'); end;
    if(Hptr >= size(H, 1) * 0.9) % if H is almost full
      %fprintf('e');
      H = [H; zeros(prealloc_size, 4)]; % add more space to H
    end
  end
  fprintf('\n');
  
  % trim invalid H (f2 should > 0, since freq bins are 1~# of freq bins)
  H = H(H(:,2)>0, : );
  
  % change the unit of f1,f2, dt, t1 to Hz and Seconds
  H = [F(H(:,1)), F(H(:,2)), H(:,3) * time_res , T(H(:,4))'];

  
  % plot peak pairs
  if DRAW 
    f1 = H(:,1); f2 = H(:,2); dt = H(:,3); t1 = H(:,4);
    plot3([t1'; (t1+dt)'], [f1'; f2'], zeros(2,length(t1)), '.-r');
    hold off;
  end
  
end


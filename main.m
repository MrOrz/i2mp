function main(i)
  %% should be calling ffmpeg to downsample here...
  wav1 = wavread(['testcases/',int2str(i),'/',int2str(i),'_ref.wav']);
  wav2 = wavread(['testcases/',int2str(i),'/',int2str(i),'_q.wav']);
  TARGET_SAMPLE_RATE = 5000;
  
  wav1 = resample(wav1(1:10*44100, 1), TARGET_SAMPLE_RATE, 44100 );
  wav2 = resample(wav2(1:10*44100, 1), TARGET_SAMPLE_RATE, 44100 );
  
  %% generating landmarks
  % lm[i,:] = ith landmark [f1, f2, dt, t1]. 
  % frequencies and times are in 'bin'.
  [lm1, F1, T1, DT1] = constellations(wav1, TARGET_SAMPLE_RATE, 1);
  fprintf('# of landmarks = %d\n',size(lm1, 1));  
  figure; % open up new figure
  [lm2, F2, T2, DT2] = constellations(wav2, TARGET_SAMPLE_RATE, 1);
  fprintf('# of landmarks = %d\n',size(lm2, 1));  
  % often the case F1 == F2, T1 == T2, DT1 == DT2.

  hashInitial;
  T = match(lm1,lm2);      
  fprintf('The time offset is %f secs.\n', T*DT2);
  

end
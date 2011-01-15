function main()
  %% should be calling ffmpeg to downsample here...
  wav = wavread('testcases/MVI_1304_ref.wav');
  
  TARGET_SAMPLE_RATE = 5000;
  
  wav = resample(wav(1:5*44100, 1), TARGET_SAMPLE_RATE, 44100 );
  %wav = resample(wav(:, 1), TARGET_SAMPLE_RATE, 44100 );
  
  %% generating landmarks
  % lm[i,:] = ith landmark [f1, f2, dt, t1]. 
  % frequencies are in Hz, times are in second.
  [lm, F, T, DT] = constellations(wav, TARGET_SAMPLE_RATE, 1);
  display(F);
  display(T);
  display(DT);

end
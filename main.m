function main()
  %% should be calling ffmpeg to downsample here...
  wav = wavread('testcases/STE-000.wav');
  
  TARGET_SAMPLE_RATE = 5000;
  
  wav = resample(wav(1:100000, 1), TARGET_SAMPLE_RATE, 44100 );
  
  %% generating landmarks
  % lm[i,:] = ith landmark [f1, f2, dt, t1]. 
  % frequencies are in Hz, times are in second.
  lm = constellations(wav, TARGET_SAMPLE_RATE, 1);
  

end
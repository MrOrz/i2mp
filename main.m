function main()
  %% should be calling ffmpeg to downsample here...
  wav1 = wavread('MVI_1304.wav');
  wav2 = wavread('MVI_1304_ref.wav');
  TARGET_SAMPLE_RATE = 5000;
  
  wav1 = resample(wav1(1:100000, 1), TARGET_SAMPLE_RATE, 44100 );
  wav2 = resample(wav2(1:100000, 1), TARGET_SAMPLE_RATE, 44100 );
  %% generating landmarks
  % lm[i,:] = ith landmark [f1, f2, dt, t1]. 
  % frequencies are in Hz, times are in second.
  lm1 = constellations(wav1, TARGET_SAMPLE_RATE, 0);
  lm2 = constellations(wav2, TARGET_SAMPLE_RATE, 1);
  hashInitial;
  match(lm1,lm2);

end
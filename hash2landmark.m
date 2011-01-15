function lm = hash2landmark(hash)
  % Input : [t,key]
  % Output : [f1,f2,dt,t]
  % 20 bits Hash value  : 8 bits of f1, 6 bits of df, 6 bits of dt

  t = hash(:,1);
  key = double(hash(:,2));
  f1 = floor(key/(2^12));
  key = key - (2^12)*f1;
  f1 = f1 + 1;  
  df = floor(key/(2^6));
  key = key - (2^6)*df;
  if df > 2^5
    df = df-2^6;
  end
  f2 = f1+df;
  dt = key;
  lm = [f1,f2,dt,t];



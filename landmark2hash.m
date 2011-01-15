function hash = landmark2hash(lm)
  % Input : [f1,f2,dt,t]
  % Output : [t,key]
  % 20 bits Hash value  : 8 bits of f1, 6 bits of df, 6 bits of dt

  t = uint32(lm(:,4));
  f1 = rem(round(lm(:,1)-1),2^8);
  df = round(lm(:,2)-lm(:,1));
  if df < 0
    df = df + 2^8;
  end
  df = rem(df,2^6);
  dt = rem(abs(round(lm(:,3))), 2^6);
  key = uint32(f1*(2^12) + df*(2^6)+dt);
  hash = [t,key];
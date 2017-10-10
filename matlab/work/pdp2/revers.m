function res = revers(n, base)


  
  res = 0.0;
  inum = n;
  power = 1.0;
  
  while (true) 
    iquot = fix(inum/base);
    irem = inum - base*iquot;
    power = power/base;
    res  =  res + irem*power;
    inum = iquot;
    if (inum < 1) return;
    end
  end

 
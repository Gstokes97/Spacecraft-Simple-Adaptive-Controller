TcSAC.Data
tauBW.Data
BW_torqueSumXYZ = sum(TcBW);
SAC_torqueSumXYZ = sum(TcSAC);

totalBW = sum(BW_torqueSumXYZ);
totalSAC = sum(SAC_torqueSumXYZ);

performanceDiff = totalBW-totalSAC

s1 = 'there is a performance difference of';
s2 = string(performanceDiff);
s3 = 'Nm by using SAC'; 

s = strcat(s1,': ',s2,s3)


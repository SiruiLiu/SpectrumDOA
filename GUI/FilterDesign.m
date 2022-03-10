Fs = 9.6e6;
FreqDiff = Fs*0.5-4.5e6;
TransitionZone = 0.07e6;
Overmeasure = 0.3e6;
Ast1 = 80;
Ap = 1;

SigFreq = Fs*0.5-FreqDiff;
Fst = (SigFreq-TransitionZone)/(Fs*0.5);


Lowpass = fdesign.highpass('N,Fp,Ast,Ap', 64, Fst, Ast1, Ap);
lpFilter = design(Lowpass, 'equiripple','SystemObject',true);
fvtool(lpFilter);
coe = lpFilter.Numerator;

save('.\FilterCOE\Coe.mat', 'coe');

function [xD]=LP_decimation_function(x,Fs,newFs)


[P,Q]=rat(newFs/Fs);

%Dise�o de un filtro �ptimo
Fsb = (Fs/Q)/2; Fpb = Fsb-500;
Rp = 1;
As = 60;
dp = (10^(Rp/20)-1)/(10^(Rp/20)+1);%par�metro de desviasi�n (banda pasante).
ds = (10^(-As/20));%par�metro de desviasi�n (banda de rechazo).
F = [Fpb Fsb]; %vector de frecuiencia (banda pasante y de rechazo).
A = [1 0]; %par�metro de amplitudes deseadas en (pb % sb).
DEV = [dp ds]; %vector de par�metros de desviaci�n.

[N,Fo,Ao,W] = firpmord(F,A,DEV,Fs); %N=190
filtroB = firpm(N,Fo,Ao,W);
xF=filter(filtroB,1,x);
xD=resample(xF,P,Q);

function [iniFin,newFs,NF]=segmentacion_Fases_Funcion(x,Fs,targetSil)

%Description:

%Main idea: indexes of the instantenous phase array are identical to the
%envelopeFiltered ones.

%

%% STEP 1: Resample at 16KHz as korean paper, by means of "LP_decimation_function".---------------------------------------------------------------------------------------------------------%       

newFs=1000;
[xD]=LP_decimation_function(x,Fs,newFs);


%% STEP 2: Calculation of hilbert transform, and envelope of new decimated signal.
xH = hilbert(xD);
xIm=imag(xH);  %Resultado Transformada de hilbert. Calculo parte imaginaria.
envelope=sqrt(xD.^2+xIm.^2); %Envelope calculation

%% STEP 3: Designing a filter which is going to be used for the segmentation process. And filtering envelope.

%Designing the segmentation filter.
filtered_band=[4 10];
[NF,b]=segmentation_Filter(newFs,filtered_band(1),filtered_band(2));

%Filtering envelope.
filteredEnvelope=filter(b,1,envelope);

%% STEP 4: Obtaining instantenous phase from the filtered envelope.

fE = hilbert(filteredEnvelope);
xImFb=imag(fE);

%Four-Quadrant Inverse Tangent. To obtain phases from the 4 quadrants
fasesFirstBand=atan2(xImFb,filteredEnvelope);    %from [-pi,pi]

%% STEP 5: Obtaining an array that contains the filtered envelope indexes of a certain phase quadrant.

quad1=fasesFirstBand>=-pi & fasesFirstBand<=-(pi/2);
quad2=fasesFirstBand>=-(pi/2) & fasesFirstBand<=-0;
quad3=fasesFirstBand>=0 & fasesFirstBand<=(pi/2);
quad4=fasesFirstBand>=(pi/2) & fasesFirstBand<=pi;

segmentedsignal=quad1+2*quad2+3*quad3+4*quad4;

% Median filtering to avoid very small transitions
segmentedsignal=medfilt1(segmentedsignal,40);
segmentedsignal=round(segmentedsignal,0);
figure; plot(segmentedsignal);

i1=[];cnt1=1;
i2=[];cnt2=1;
i3=[];cnt3=1;
i4=[];cnt4=1;

for i=1:length(segmentedsignal)
    if segmentedsignal(i)==1
        i1(cnt1)=i;
        cnt1=cnt1+1;
    end
    if segmentedsignal(i)==2
        i2(cnt2)=i;
        cnt2=cnt2+1;
    end
    if segmentedsignal(i)==3
        i3(cnt3)=i;
        cnt3=cnt3+1;
    end
    if segmentedsignal(i)==4
        i4(cnt4)=i;
        cnt4=cnt4+1;
    end 
end

%% Step 6: Get the begining and end point of each segment.

transiciones=find(diff(segmentedsignal)~=0);
nSeg=length(transiciones)+1;

inicio=zeros(nSeg,1);
final=zeros(nSeg,1);

for i=1:nSeg
    if i>1 && i~=nSeg
        final(i)=transiciones(i);
        inicio(i)=final(i-1)+1;
    end
    if i==1
        final(i)=transiciones(i);
        inicio(i)=1;
    end
    if i==nSeg
        final(i)=length(xD);
        inicio(i)=final(i-1)+1;
    end
end

iniFin=[inicio final]; %¿Pos-procesado segmentos sin contenido?

%% STEP 7: Several Displays


%Display of filtered envelope segments:
% filt_Env_Segm_Plot(i1,i2,i3,i4,newFs,filteredEnvelope);

%Display1 of original signal segments:
%x_Segm_Plot1(i1,i2,i3,i4,newFs,Fs,x);

%Display2 of original signal segments:
%Visualizacion: Señal original, barra de segmentos, lineas verticales
%delimitadoras.
x_Segm_Plot2(i1,i2,i3,i4,x,newFs,Fs,targetSil,iniFin);

%% STEP 8: Functions

function [NF,b]=segmentation_Filter(newFs,Fp1,Fp2)
%[z,p,k]=ellip(3,5,40,bandaTheta*2/Fs);
%[b,a]=zp2tf(z,p,k);
Fst1 = Fp1-3;
Fst2 = Fp2+3;
Ap=6;
Ast1=40;
Ast2=40;
Hf = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2,newFs);
H2= design(Hf,'equiripple');
NF=length(H2.numerator);
b=H2.numerator;
end

function [i1r, i2r, i3r, i4r]=samples_Del_newFs2Fs(x,i1, i2, i3, i4, newFs, Fs)

%With the real Fs signal:
i1r=i1*(Fs/newFs);
if i1r(end)>length(x)
    i1r(end)=length(x);
end

i2r=i2*(Fs/newFs);
if i2r(end)>length(x)
    i2r(end)=length(x);
end

i3r=i3*(Fs/newFs);
if i3r(end)>length(x)
    i3r(end)=length(x);
end

i4r=i4*(Fs/newFs);
if i4r(end)>length(x)
    i4r(end)=length(x);
end
end


function filt_Env_Segm_Plot(xD,i1,i2,i3,i4,newFs,filteredEnvelope)

tr=(1:length(xD))/newFs;
figure; plot(tr,filteredEnvelope); hold on;
if isempty(i1)
else
    plot(i1/newFs,filteredEnvelope(i1),'.','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5); hold on;
end
if isempty(i2)
else
    plot(i2/newFs,filteredEnvelope(i2),'.','MarkerEdgeColor','b','MarkerFaceColor','r','MarkerSize',5); hold on;
end
if isempty(i3)
else
    plot(i3/newFs,filteredEnvelope(i3),'.','MarkerEdgeColor','y','MarkerFaceColor','r','MarkerSize',5); hold on;
end
if isempty(i4)
else
    plot(i4/newFs,filteredEnvelope(i4),'.','MarkerEdgeColor','g','MarkerFaceColor','r','MarkerSize',5); hold off;
end

end

function x_Segm_Plot1(i1,i2,i3,i4,newFs,Fs,x)

[i1r, i2r, i3r, i4r]=samples_Del_newFs2Fs(x,i1, i2, i3, i4, newFs, Fs);

tx=(1:length(x))./Fs;
figure; plot(tx,x,'k'); hold on;
if isempty(i1r)
else
    plot(i1r/Fs,x(i1r),'.','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',15); hold on;
end
if isempty(i2r)
else
    plot(i2r/Fs,x(i2r),'.','MarkerEdgeColor','b','MarkerFaceColor','r','MarkerSize',15); hold on;
end
if isempty(i3r)
else
    plot(i3r/Fs,x(i3r),'.','MarkerEdgeColor','y','MarkerFaceColor','r','MarkerSize',15); hold on;
end
if isempty(i4r)
else
    plot(i4r/Fs,x(i4r),'.','MarkerEdgeColor','g','MarkerFaceColor','r','MarkerSize',15); hold off;
end
    
end

function x_Segm_Plot2(i1,i2,i3,i4,x,newFs,Fs,targetSil,iniFin) %Fs=16000, newFs=1000

[i1r, i2r, i3r, i4r]=samples_Del_newFs2Fs(x,i1, i2, i3, i4, newFs, Fs);

figure;
minValueX=min(x);
minValueX1=repelem(minValueX,length(i1r));
minValueX2=repelem(minValueX,length(i2r));
minValueX3=repelem(minValueX,length(i3r));
minValueX4=repelem(minValueX,length(i4r));
tx=(1:length(x))./Fs;
plot(tx,x,'k');hold on;
plot(i1r/Fs,minValueX1,'s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10);hold on;
plot(i2r/Fs,minValueX2,'s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10);hold on;
plot(i3r/Fs,minValueX3,'s','MarkerEdgeColor','y','MarkerFaceColor','y','MarkerSize',10);hold on;
plot(i4r/Fs,minValueX4,'s','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',10);hold on;
xlabel('Tiempo(s)'); ylabel('Amplitud(Pa)'); title(['Silaba "',targetSil,'" señal original sin decimar'],'Interpreter' ,'none');

%lineas verticales delimitadoras
del=iniFin(:,2); %array con todos los finales de segmento.
del=del*(Fs/newFs); %Esas muestras les cambio su Fs por la original
ax=gca;
ylims=ax.YLim;  %Te da un array de los limites en valor del plot
ax.YLim=[(minValueX+minValueX*0.8) ylims(2)];
for g=1:length(del)
    x= repelem(del(g)/Fs,length(minValueX:0.01:ylims(2)));
    plot(x,minValueX:0.01:ylims(2),'b','LineStyle','--','LineWidth',2); hold on;
end

end


end
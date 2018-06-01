function [x_framed,splitter,iniFinS]=framer(x,Fs,L,R,wt)

%x_framed: tiene filas con zeros al final posiblemente, a la hora de
%seleccionar puramente los valores de un frame: x_framed(i,x_framed(1,x_framed(1,:)~=0);
%O tambien: x_framed(1,1:L);

M=length(x);
t=(1:1:M)/Fs;

w1=window(str2func(wt),L);   %Crea ventana en columna                             


%MODIFICAR EST DE ABAJO
if mod(M,R)==0 && R~=0
    splitter=1:R:M-L+1;
    inii=splitter';
    fini=splitter'+L-1;
    iniFinS=[inii fini];
    nSeg=length(iniFinS(:,1));
    x_framed=zeros(nSeg,L-1);
    [m,n]=size(iniFinS);
    iniFinS(m,n)=length(x);
    for i=1:1:nSeg
        x_framed(i,1:L)=w1.*x(splitter(i):splitter(i)+L-1);
    end
else
    if R~=0
        splitter=1:R:M-L+1;
        inii=splitter';
        fini=splitter'+L-1;
        iniFinS=[inii fini];
        nSeg=length(iniFinS(:,1));
        x_framed=zeros(nSeg,L-1);
        [m,n]=size(iniFinS);
        iniFinS(m,n)=length(x);
        for i=1:1:nSeg-1 
            x_framed(i,1:L)=w1.*x(iniFinS(i,1):iniFinS(i,2)); 
        end
         x_framed(i+1,1:length(x(iniFinS(i+1,1):iniFinS(i+1,2))))=x(iniFinS(i+1,1):iniFinS(i+1,2));
    end
end



if mod(M,R)==0 && R==0  %Hay una segmentacion sin solapacion, y es par
    splitter=1:L-1:M;
    splitter=[splitter M];
    for i=1:1:length(splitter)-1
        x_framed(i,1:length(x(splitter(i):splitter(i+1))))=w1.*x(splitter(i):splitter(i+1));
    end
   
elseif R==0 %Sin solapacion y no es par
     
        splitter=1:L-1:M;   %Clave L-1
        splitter=[splitter M];
        nSeg=length(splitter)-1;
        iniFinS=zeros(nSeg,2);
        for k=1:nSeg
            if k==1
                iniFinS(k,1:2)=[splitter(1,k) splitter(1,k+1)]; %when k=1
            end
            if k>1 && k<nSeg
                iniFinS(k,1:2)=[(iniFinS(k-1,2)+1) (splitter(1,k+1)+k-1)];
            end
            if k==nSeg
                iniFinS(k,1:2)=[(iniFinS(k-1,2)+1) splitter(1,k+1)]; %when k==NSeg
            end
        end
     
        x_framed=zeros(nSeg,L-1);
        
        for i=1:1:nSeg-1
            x_framed(i,1:L)=x(iniFinS(i,1):iniFinS(i,2));
        end
        x_framed(nSeg,1:length(x(iniFinS(nSeg,1):iniFinS(nSeg,2))))=x(iniFinS(nSeg,1):iniFinS(nSeg,2));
       
end



    
    
    
    
    
    
    
    
    
    
    
    
    


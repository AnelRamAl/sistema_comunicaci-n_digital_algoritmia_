clear all
clc

%MUESTREO
Fs=8e3; 
n_bits=16; %bits por muestra
canal=1; 
tiempo=2; %duración de audio

recObj = audiorecorder(Fs,n_bits,canal);
disp('Start speaking.')
recordblocking(recObj,tiempo);
disp('End of Recording.');
myRecording = getaudiodata(recObj);
save('audio_de_prueba.mat','myRecording');
%recordblocking(recObj,Fs)%control de tiempo hasta q acabe el audio
load audio_de_prueba.mat %generación del archivo de audio wav
filename = 'audio_de_prueba.wav'; 
y=myRecording;
audiowrite(filename,y,Fs);
clear y

%Señal muestreada
[y,Fm] = audioread('audio_de_prueba.wav'); 
sound(y,Fm)
figure(1)
subplot(1,2,1);plot(y);ylabel('F(t)[volts]');xlabel('t(seg)');title('a)Señal muestreada (t)')
[b1,a1]=hist(y,100);
s0=sort(y);
u0=unique(s0,'stable');
z0=b1/sum(b1);
subplot(1,2,2);bar(a1,z0);ylabel('Fx(xi)');xlabel('x');title('b)PDF de la señal muestreada')
fun_A=0;
fun_miu=0;
parametro=input('Ingrese opción de compresión ley A (ESTÁNDAR)=1, Ley B (ESTANDAR)=2, Ley A=3 , Ley B=4:     ')
switch parametro
    case 1%ley A
        A = 87.6;
        fun_A=1;
    case 2%ley B
        miu=255;
        fun_miu=1;
    case 3%ley A
        A=input('Ingrese parametro A')
        fun_A=1;
    otherwise 
        miu=input('Ingrese parametro miu')
        fun_miu=1;
end
 
%COMPANDING ley A
if fun_A==1;
%A = 87.6;
v=max(y); 
%comp_y = compand(y,A,v,'a/compressor');
total=tiempo*Fs; 
ma=max(abs(y));
ganancia=1/ma;
y=y*ganancia;
i=1:1:total;
 comp_y=sign(y(i)).*((((abs(y(i))<=(v/A))&(abs(y(i))>=0)).*((A*abs(y(i))/(1+log(A)))))+(((abs(y(i))<=v)&(abs(y(i))>(v/A))).*((v*(1+log(A*abs(y(i))/v))/(1+log(A))))));   
%sound(comp_y,Fm)
comp_x=hist(comp_y,100);
%COMPANDING ley miu
figure(2)
subplot(1,2,1);plot(comp_y);ylabel('F(t)');xlabel('t');title('c)Ley A/miu señal en el tiempo')
subplot(1,2,2);bar(comp_x);ylabel('Fx(xi)');xlabel('x');title('d)Ley A/miu PDF')

elseif fun_miu==1;
%miu=255;
v=max(y);   
ma=max(abs(y));
total=tiempo*Fs; 
ganancia=1/ma;
y=y*ganancia;%%%%COMENTAAAR
i=1:1:total;
    comp_y = sign(y(i)).*(v*log(1+miu*abs(y(i))/v)/log(1+miu));            %comp_u_y = compand(y,miu,v,'mu/compressor')
%sound(comp_u_y,Fm)
comp_u_x=hist(comp_y,100);
figure(2)
subplot(1,2,1);plot(comp_y);ylabel('F(t)');xlabel('t');title('c)Ley A/miu señal en el tiempo')
subplot(1,2,2);bar(comp_u_x);ylabel('Fx(xi)');xlabel('x');title('d)Ley A/miu PDF')

else error=0
end

%CUANTIZACIÓN UNIFORME
niveles=256; 
mem=0;
cuant=0;
q=(max(comp_y)-min(comp_y))/niveles; 
ganancia=1/((max(comp_y)+(2*niveles-1)*(q/2))-(min(comp_y)+(q/2)));
gan=zeros(1,total);
for i=1:1:total;
    
    for n=1:1:niveles;    
    q_n= min(comp_y)+(2*n-1)*(q/2);
     int_inf=min(comp_y)+(n-1)*q;
     int_sup=min(comp_y)+(n)*q;
       cuant=((comp_y(i)>=int_inf)&(comp_y(i)<int_sup)).*q_n;
       mem=mem+cuant;
       if cuant~=0;
           n=niveles;
           cuant=0;
       end
    end 
gan(i)=mem; 
mem=0;
end
figure(3)
plot(gan);ylabel('F(t)');xlabel('t');title('f)Cuantización')

%CODIFICACIÓN
bits=round(log2(niveles));
larg=2^bits; 
mem=0;
X=0;
h=zeros(1,total);
for i=1:1:total;
   for n=1:1:niveles;   
    q_n= min(comp_y)+(2*n-1)*(q/2);    
      if gan(i)==q_n;
        mem=dec2bin(n-1); %vector de caracteres
        X=str2num(mem);
        n=niveles;
      end
   end
 h(i)=X;
end
% figure(4)
% plot(h)
% set(gca,'ytick', [000 001 010 011 100 101 110 111]); 
% ylabel('F(t)');xlabel('t');title('G)Codificación')


cadena=zeros(1,bits*total);%lenght(h)=total=frec*tiempo
posicion=1;
for i=1:1:total;
    mod_n=h(i);
    for m=bits:-1:2;%-1 posición de la asignación del modulo
        div=floor(mod_n/(10^(m-1)));%cociente
        mod_n=mod(mod_n,10.^(m-1));
        
                if 10^(m-1)==10;% modulo
                    cadena(posicion)=div;
                    cadena(posicion+1)=mod_n;
                    posicion=posicion+1;
                else cadena(posicion)=div;
                end 
         posicion=posicion+1;
    end
end
% figure(5)
% plot(cadena)

%CODIGO DE LINEA
codigo=input('Seleccione código de linea:  1=Manchester    2=NRZ-L:          ')          
switch codigo
    case 1
%cadena=[1 0 0 0 0]
T=4;    %T=periodo  %flanco de subida 0 flanco de bajada 1 
A_menos = [-1*ones(1,T/2) ones(1,T/2)];
A_mas = [ones(1,T/2) -1*ones(1,T/2)];
lv=bits*total;%lv=length(cadena);
salida1=[];
    for a=1:lv
        if (cadena(a)== 1)
            salida1=[salida1 A_mas];
        else 
            salida1=[salida1 A_menos];
        end
    end

    case 2
        lv=bits*total;%lv=length(cadena);
        T=2;
   %lv=bits*total;
   lv=length(cadena);
        salida=[];
        unos=ones(1,T);
        for a=1:lv
            if(cadena(a)==1)
                salida=[salida unos];
            else
                salida=[salida -unos];
            end
        end    
        
    otherwise  cod=0; %codigo 3
        
end
        

%y = awgn(salida1,snr,sigpower)%The scalar snr specifies the signal-to-noise ratio per sample, in dB.
                              %accepts the input sigpower, which specifies the power of x in dBW.
snr=input('Ingrese SNR en dB:               ')          
y = awgn(salida1,snr); %This syntax assumes that the power of x=salida1 is 0 dBW.                              
figure(6)
plot(y)
hold on 
plot(salida1)
ylabel('F(t)');xlabel('t');title('4.Código de linea y codigo + ruido AGWN')

%%Receptor%%
%RECEPTOR DEL CODIGO DE LINEA
if codigo==1
 m=sign(salida1(1:4:end));
 receptor= gt(m,0);
 
elseif codigo==2
     salida=gt(salida(1:2:end),0);
end

%PRE-DECODIFICACIÓN %%pasar bits a datos 111 101 111
posicion=1;
mem=0;
y=zeros(1,total);%nueva cadena de datos para la decodificación
for i=1:1:total-1;
    for n=bits:-1:1;%-1 posición de la asignación del modulo
        numero=receptor(posicion)*10.^(n-1);
        mem=numero+mem;
        posicion=posicion+1;
    end 
    y(i)=mem;
    mem=0;
end
% figure(7)
% plot(y)

%DECODIFICACIÓN
bits=round(log2(niveles));
larg=2^bits; %num de cod     larg==niv 
mem_r=0;
q_n_r=zeros(1,total);
for i=1:1:total;
   for n=1:1:niveles;  
       mem_r=dec2bin(n-1);
       X=str2num(mem_r);
      if y(i)==X;
       buf= min(comp_y)+(2*n-1)*(q/2);  
       n=niveles;
      end
   end
 q_n_r(i)=buf;
end
figure(7)
plot(q_n_r);ylabel('F(t)');xlabel('t');title('6.a)decodificación')

%CONTROL DE GANANCIA

gan_r=zeros(1,total);
for i=1:1:total;
    gan_r(i)=q_n_r(i)*ganancia;
end
figure(8)
plot(gan_r);ylabel('F(t)');xlabel('t');title('6.b)control de ganacia')

%DECOMPANDING ley miu
if fun_miu==1;
   expan_u_y = compand(q_n_r,miu,v,'mu/expander');
i=1:1:total;
    %expan_u_y=sign(gan_r(i)).*((v/miu)*(((1+miu).^(abs(gan_r(i))/v))-1));
expan_u_x=hist(expan_u_y,100);
sound(expan_u_y,Fm)
figure(9)
subplot(1,2,1);plot(expan_u_y);ylabel('F(t)');xlabel('t');title('6.e)Señal de voz recuperada')
subplot(1,2,2);bar(expan_u_x);ylabel('Fx(xi)');xlabel('x');title('6.d)PDF de la voz recuperada')
%DECOMPANDING ley a
elseif fun_A==1;
   %expan_y = compand(in,87.6,v,'a/expander');
i=1:1:total;
   expan_y=sign(gan_r(i)).*((((abs(gan_r(i))<=(v/(1+log(A))))&(abs(gan_r(i))>=0)).*abs(gan_r(i))*((1+log(A))/A))+(((abs(gan_r(i))<=v)&(abs(gan_r(i))>=(v/(1+log(A))))).*((v*(exp(abs(gan_r(i))*((1+log(A))/v)-1))/A))));   
expan_x=hist(expan_y,100);
sound(expan_y,Fm)
figure(9)
subplot(1,2,1);plot(expan_y);ylabel('F(t)');xlabel('t');title('6.e)Señal de voz recuperada')
subplot(1,2,2);bar(expan_x);ylabel('Fx(xi)');xlabel('x');title('6.d)PDF de la voz recuperada')
else error=0
end
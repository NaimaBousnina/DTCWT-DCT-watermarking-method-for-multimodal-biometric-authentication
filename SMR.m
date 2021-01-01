function ZO_Polaire = SMR()

%loading the database of minutiaes.
load('Or_Minutias.mat'); 
for i=1:length(Minutia)
nb_Minu_image{i}=size(Minutia{i}); % Number of minutiae of each image.
end

Taille_image=[560,296];
sigmaN=3.2;
sigmaL=0.32;
sigmaO=3.87;
lambdaL=0.001;
lambdaH=0.53;

%Location-based spectral minutiae representation(SML) [domaine de Fourier]
windowSize=5;
for i=1:length(Minutia)
    Flx{i}=linspace(-windowSize,windowSize,Taille_image(2));
    Fly{i}=linspace(-windowSize,windowSize,Taille_image(1));
    [wx{i},wy{i}]=meshgrid(Flx{i},Fly{i});
    ZL{i} = zeros(Taille_image(1),Taille_image(2));
    ZL{i} = complex(ZL{i},0);
    for j=1:nb_Minu_image{i}(1)
        ZL{i}=ZL{i}+ exp(-1j*((wx{i}*Minutia{i}(j,1))+(wy{i}*Minutia{i}(j,2))));
    end
    D1=2*(sigmaL^(-2));
    D2=exp(-((wx{i}).^2+(wy{i}).^2)/D1);
    ZL{i}=ZL{i}.*D2;
     
   %plot LSMR 
   %figure(i);
   %imagesc(abs(ZL{i})); title('SML dans la domaine de Fourier');
end

%Orientation-based spectral minutiae representation(SMO)[domaine de Fourier]
for i=1:length(Minutia)
    Flx{i}=linspace(-windowSize,windowSize,Taille_image(2));
    Fly{i}=linspace(-windowSize,windowSize,Taille_image(1));
    [wx{i},wy{i}]=meshgrid(Flx{i},Fly{i});
    ZO{i} = zeros(Taille_image(1),Taille_image(2));
    ZO{i} = complex(ZO{i},0);
    for j = 1:nb_Minu_image{i}(1)
        ZO{i}= ZO{i} + 1j*(wx{i}*cos(Minutia{i}(j,4))+wy{i}*sin(Minutia{i}(j,4))) .* exp(-1j*(wx{i}*Minutia{i}(j,1)+wy{i}*Minutias{i}(j,2))) ;
    end
    D3=exp(-((wx{i}.^2)+(wy{i}.^2))/(2*(sigmaO^(-2))));
    ZO{i}=ZO{i}.*D3;
    
  %Plot OSMR
  %figure(i);
  %imagesc(abs(ZO{i})); title('SMO dans le domane de Fourier');
end

%Plot in a Log-polaire space.
%In the radial direction M = 128 samples are used between lambdaL and lambdaH
%In the angular direction N = 256 samples are used between  0 and 2pi

for i=1:length(Minutia)
    rcoords = linspace(lambdaL,lambdaH,128);
    thcoords = linspace(0,pi,256);
    [ri{i},thi{i}] = meshgrid(rcoords,thcoords);
    [Lx{i},Ly{i}] = pol2cart(thi{i},ri{i});
      
    %Plot the SML in a Log-polaire space.
    ZL_Plaire{i} = interp2(wx{i},wy{i},abs(ZL{i}),Lx{i},Ly{i});
    ZL_Plaire{i}=ZL_Plaire{i}';
    %figure(i);
    %imagesc(ZL_Plaire{i}); title('Les coordonnees polaire de SML');
      
    %Plot the SMO in a Log-polaire space.
    ZO_Polaire{i} = interp2(wx{i},wy{i},abs(ZO{i}),Lx{i},Ly{i});
    ZO_Polaire{i}=ZO_Polaire{i}';
    %figure(i);
    %imagesc(ZO_Polairee{i}); title('Les coordonnees polaires de SMO');
end

end

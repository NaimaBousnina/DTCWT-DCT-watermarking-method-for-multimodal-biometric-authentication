clear all;
close all;
clc;

     %%%% Watermark preparation %%%%
     
% Getting the fingerprint features 
ZO_Polaire = SMR() %The spectral minutiae representation
Vector_Reduced = CPCA(ZO_Polaire) %Apply the Column Principal Component Analysis (CPCA) feature 
                                  %reduction method to decorrelate features and concentrate power.
% Normlize the values of the reduced template (for the quantization).
for k=1:length(Vector_Reduced)
    moyenne{k}=mean(Vector_Reduced{k});
    Norm_Template{k}=Vector_Reduced{k}/moyenne{k};
end
    Binary_Minutia_all = Quantise(Norm_Template) %Quantize the reduced features to get the watermark.

    %%%% Watermark embedding %%%%
   
load('Face_im.mat');
for k=1:length(faces)
    %Apply first level DTCWT to decompose the face image
    x{k}=double(faces{k});
    [Yl{k},Yh{k}] = dtwavexfm2(x{k},1,'near_sym_a','qshift_a');  

    %Take Two high frequency sub_bands.
     One_high_freq{k}=Yh{k}{1}(:,:,1);
     Two_high_freq{k}=Yh{k}{1}(:,:,4);
       
    %devide the sub bands frequency into 4*4 blocs.
     N_L{k} = 4*ones(1,64);
     N_C{k} = 4*ones(1,64);
     One_High_freq_blocs{k} =mat2cell(One_high_freq{k},N_L{k},N_C{k});
     Two_High_freq_blocs{k} =mat2cell(Two_high_freq{k},N_L{k},N_C{k});
       
     %Apply DCT to each block of sub bands.
     for m=1:size(One_High_freq_blocs{k},1)
         for n=1:size(One_High_freq_blocs{k},2)
             One_High_freq_blocs_dct{k}{m,n}=dct2(One_High_freq_blocs{k}{m,n});
             Tow_High_freq_blocs_dct{k}{m,n}=dct2(Two_High_freq_blocs{k}{m,n});
         end
     end
end

[PN_1 PN_0 Z]=water_embed(One_High_freq_blocs_dct,Tow_High_freq_blocs_dct,Binary_Minutia_all,Yh,Yl);

      %%%% Watermark extraction %%%%  
Extracted_water = water_Extraction(Z,Binary_Minutia_all,PN_0,PN_1);

%%%% END %%%%

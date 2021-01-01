function Extracted_watermark_bits = water_Extraction(Z,Binary_Minutia_all,PN_0,PN_1)

midband=[0,0,1,1;   
         0,1,1,0;
         1,1,0,0;
         1,0,0,0];
blocksize=4;
for k=1:length(Z)
    Extracted_watermark_bits{k}=[];
    %Apply first level DTCWT to the watermarked image.
    [Yl_water_im{k},Yh_water_im{k}] = dtwavexfm2(Z{k},1,'near_sym_a','qshift_a');
        
    %Take two high frequency sub_bands.
    One_high_water_im{k}=Yh_water_im{k}{1}(:,:,1);
    Two_high_water_im{k}=Yh_water_im{k}{1}(:,:,4);
           
    %Devide the sub band frequencies into 4*4 blocs.
    N_L = 4*ones(1,64);
    N_C = 4*ones(1,64);
    One_High_blocs_water_im{k} =mat2cell(One_high_water_im{k},N_L,N_C);
    Two_High_blocs_water_im{k} =mat2cell(Two_high_water_im{k},N_L,N_C);
           
    %Apply DCT to each block of sub bands.
    for m=1:size(One_High_blocs_water_im{k},1)
        for n=1:size(One_High_blocs_water_im{k},2)
            One_High_blocs_water_im_dct{k}{m,n}=dct2(One_High_blocs_water_im{k}{m,n});
            Two_High_blocs_water_im_dct{k}{m,n}=dct2(Two_High_blocs_water_im{k}{m,n});
        end
    end
         
    %Extract the mid band coefficients
    for m=1:size(One_High_blocs_water_im{k},1)
        for n=1:size(One_High_blocs_water_im{k},2)
            for r=1:blocksize
                for l=1:blocksize
                    if (midband(r,l)==1)
                        One_High_mid_band_coef{k}{m,n}(r,l)=One_High_blocs_water_im_dct{k}{m,n}(r,l);
                        Two_High_mid_band_coef{k}{m,n}(r,l)=Two_High_blocs_water_im_dct{k}{m,n}(r,l);
                    end
                 end
            end
        end
     end
         
    for m=1:size(One_High_blocs_water_im{k},1)
        for n=1:size(One_High_blocs_water_im{k},2)
            One_High_mid_band_coef_prim{k}{m,n}=One_High_mid_band_coef{k}{m,n}';
            One_High_mid_band_coef_vec{k}{m,n} = One_High_mid_band_coef_prim{k}{m,n}(One_High_mid_band_coef_prim{k}{m,n}~=0);
                 
            Two_High_mid_band_coef_prim{k}{m,n}=Two_High_mid_band_coef{k}{m,n}';
            Two_High_mid_band_coef_vec{k}{m,n} = Two_High_mid_band_coef_prim{k}{m,n}(Two_High_mid_band_coef_prim{k}{m,n}~=0);
        end
    end
         
    %Correlation computation.
    for m=1:size(One_High_blocs_water_im{k},1)
        for n=1:size(One_High_blocs_water_im{k},2)  
            One_Corr_High_0{k}{m,n}=corr(real(One_High_mid_band_coef_vec{k}{m,n}), PN_0{k}');
            One_Corr_High_1{k}{m,n}=corr(real(One_High_mid_band_coef_vec{k}{m,n}), PN_1{k}');
                 
            Two_Corr_High_0{k}{m,n}=corr(real(Two_High_mid_band_coef_vec{k}{m,n}), PN_0{k}');
            Two_Corr_High_1{k}{m,n}=corr(real(Two_High_mid_band_coef_vec{k}{m,n}), PN_1{k}');
         end
    end
         
    %Extract the watermark bits.
    One_indice_wat_bit{k}=1;
    Two_indice_wat_bit{k}=1;
        
    for m=1:size(One_High_blocs_water_im{k},1)
        for n=1:size(One_High_blocs_water_im{k},2)
            
            %First high sub bands.
            if One_Corr_High_0{k}{m,n}>=One_Corr_High_1{k}{m,n}     
               One_Extracted_watermark_bits{k}(One_indice_wat_bit{k})=0;
               One_indice_wat_bit{k}=One_indice_wat_bit{k}+1;
            else
               One_Extracted_watermark_bits{k}(One_indice_wat_bit{k})=1;
               One_indice_wat_bit{k}=One_indice_wat_bit{k}+1;
            end   
            
            %Second high sub bands.
            if Two_Corr_High_0{k}{m,n}>=Two_Corr_High_1{k}{m,n}     
               Two_Extracted_watermark_bits{k}(Two_indice_wat_bit{k})=0;
               Two_indice_wat_bit{k}=Two_indice_wat_bit{k}+1;
            else
               Two_Extracted_watermark_bits{k}(Two_indice_wat_bit{k})=1;
               Two_indice_wat_bit{k}=Two_indice_wat_bit{k}+1;
            end
         end
    end
         
    Extracted_watermark_bits{k}=[One_Extracted_watermark_bits{k} Two_Extracted_watermark_bits{k}];
end

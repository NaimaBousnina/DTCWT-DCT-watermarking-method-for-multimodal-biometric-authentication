function [PN_1 PN_0 Z]=water_embed(One_High_freq_blocs_dct,Tow_High_freq_blocs_dct,Binary_Minutia_all,Yh,Yl)

%Define the mid-band frequencies of an 4x4 dct
midband=[0,0,1,1;   
         0,1,1,0;
         1,1,0,0;
         1,0,0,0];
blocksize=4;
for k=1:length(One_High_freq_blocs_dct)    
    %Generate the two pseudo-random sequences PN0 and PN1.
    PN_1{k} = rand(1,sum(sum(midband)));
    ind{k} = PN_1{k} >= 0.5;
    PN_1{k}(ind{k}) = 1;
    PN_1{k}(~ind{k}) = -1;
    
    PN_0{k} = rand(1,sum(sum(midband)));
    ind{k} = PN_0{k} >= 0.5;
    PN_0{k}(ind{k}) = 1;
    PN_0{k}(~ind{k}) = -1;
                    
    %Devide the watermark into two blocks.
    One_Binary_Minutia_all{k}=Binary_Minutia_all{k}(1:4096);
    Two_Binary_Minutia_all{k}=Binary_Minutia_all{k}(4097:8192);
                
                %%%%%%+++++++++++++++++++++++++++++++++++++++++%%%%%%%%
                %%%%%% Embed into the first high DCT sub bands %%%%%%%%
                %%%%%%+++++++++++++++++++++++++++++++++++++++++%%%%%%%%
                
    ind_water_bit{k}=1; % Watermark bits' indix.
    for m=1:size(One_High_freq_blocs_dct{k},1)
        for n=1:size(One_High_freq_blocs_dct{k},2)
            if (ind_water_bit{k}<=length(One_Binary_Minutia_all{k}))&(One_Binary_Minutia_all{k}(ind_water_bit{k})==0) % if message bit contains zero then embed pn_sequence_zero into the mid-band componants.
               ll=1; %Index of the seuence PN0 elements.
               one_High_SB{k}{m,n}=One_High_freq_blocs_dct{k}{m,n};
               for r=1:blocksize
                   for l=1:blocksize
                       if (midband(r,l)==1)
                           one_High_SB{k}{m,n}(r,l)=One_High_freq_blocs_dct{k}{m,n}(r,l)+20*PN_0{k}(ll);
                           ll=ll+1;
                       end
                   end
               end
               ind_water_bit{k}=ind_water_bit{k}+1;
             %Otherwise, embed pn_sequence_one into the mid-band componants.  
             elseif (ind_water_bit{k}<=length(One_Binary_Minutia_all{k}))&(One_Binary_Minutia_all{k}(ind_water_bit{k})==1)
                    ll=1; %Index of the seuence PN1 elements.
                    one_High_SB{k}{m,n}=One_High_freq_blocs_dct{k}{m,n};
                    for r=1:blocksize
                        for l=1:blocksize
                            if (midband(r,l)==1)
                               one_High_SB{k}{m,n}(r,l)=One_High_freq_blocs_dct{k}{m,n}(r,l)+20*PN_1{k}(ll);
                               ll=ll+1;
                            end
                        end
                    end
                    ind_water_bit{k}=ind_water_bit{k}+1;
             else
                    one_High_SB{k}{m,n}=One_High_freq_blocs_dct{k}{m,n};
            end     
         end
   end 
                       
                %%%%%%+++++++++++++++++++++++++++++++++++++++++%%%%%%%%
                %%%%%% Embed into the second high DCT sub bands %%%%%%%%
                %%%%%%+++++++++++++++++++++++++++++++++++++++++%%%%%%%%
                
       ind_water_bit{k}=1; %Watermark bits' indix.
       for m=1:size(Tow_High_freq_blocs_dct{k},1)
           for n=1:size(Tow_High_freq_blocs_dct{k},2)
               if (ind_water_bit{k}<=length(Two_Binary_Minutia_all{k}))&(Two_Binary_Minutia_all{k}(ind_water_bit{k})==0) % if message bit contains zero then embed pn_sequence_zero into the mid-band componants.
                  ll=1; % %Index of the seuence PN0 elements.
                  Tow_High_SB{k}{m,n}=Tow_High_freq_blocs_dct{k}{m,n};
                  for r=1:blocksize
                      for l=1:blocksize
                          if (midband(r,l)==1)
                             Tow_High_SB{k}{m,n}(r,l)=Tow_High_freq_blocs_dct{k}{m,n}(r,l)+30*PN_0{k}(ll);
                             ll=ll+1;
                          end
                      end
                   end
                   ind_water_bit{k}=ind_water_bit{k}+1;
               %Otherwise, embed pn_sequence_one into the mid-band componants.
               elseif (ind_water_bit{k}<=length(Two_Binary_Minutia_all{k}))&(Two_Binary_Minutia_all{k}(ind_water_bit{k})==1)                           
                      ll=1; %Index of the seuence PN1 elements.
                      Tow_High_SB{k}{m,n}=Tow_High_freq_blocs_dct{k}{m,n};
                      for r=1:blocksize
                          for l=1:blocksize
                              if (midband(r,l)==1)
                                 Tow_High_SB{k}{m,n}(r,l)=Tow_High_freq_blocs_dct{k}{m,n}(r,l)+30*PN_1{k}(ll);
                                 ll=ll+1;
                              end
                          end
                      end
                      ind_water_bit{k}=ind_water_bit{k}+1;
               else
               Tow_High_SB{k}{m,n}=Tow_High_freq_blocs_dct{k}{m,n};
               end                              
            end
      end                 
end

for k=1:length(One_High_freq_blocs_dct)
    %Apply IDCT (inverse DCT) to each high sub band block
    for m=1:size(one_High_SB{k},1)
        for n=1:size(one_High_SB{k},2)
            One_High_SB_IDCT{k}{m,n}=idct2(one_High_SB{k}{m,n}); %%% First high subband.
            Tow_High_SB_IDCT{k}{m,n}=idct2(Tow_High_SB{k}{m,n}); %%% Second high subband.
        end
     end
        
     %Cell array to ordinary array.
     First_High_SB{k}=cell2mat(One_High_SB_IDCT{k});
     Second_High_SB{k}=cell2mat(Tow_High_SB_IDCT{k});
        
     %Embed the High sub bands into the global image
     Yh_After{k}=Yh{k};
     Yh_After{k}{1}(:,:,1)=First_High_SB{k};
     Yh_After{k}{1}(:,:,4)=Second_High_SB{k};
        
     %Apply the inverse DTCWT (IDTCWT)
     Z{k} = dtwaveifm2(Yl{k},Yh_After{k},'near_sym_a','qshift_a');
     %imwrite(uint8(Z{k}),['Water_Image\IM' num2str(k) '.tif'],'tif');
end

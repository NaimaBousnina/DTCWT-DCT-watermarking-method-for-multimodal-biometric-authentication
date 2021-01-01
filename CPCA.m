function Vector_Reduced = CPCA(ZO_Polaire)

nb_representations=length(ZO_Polaire); 
k=1;
m=1;
nb_im_per_user=8;  %Number of images per user.
ZO{m}=zeros(256,128);
while k<=nb_representations
      if k<=nb_im_per_user
         ZO_Pol{k}=ZO_Polaire{k}';
         ZO{m}=[ZO{m};ZO_Pol{k}];
         k=k+1;
      else
         nb_im_per_user=nb_im_per_user+8;  
         m=m+1;
         ZO{m}=zeros(256,128);
      end
end

nb_im_per_user=8;
for i=1:(nb_representations/nb_im_per_user)
    ZO{i}=ZO{i}';
    ZO{i}=ZO{i}(:,257:size(ZO{i},2));
    
    %Compute the mean of eac representation.
    moyenne1{i}=mean(ZO{i});
    moyenne2{i}=mean(moyenne1{i});
    
    %sustract the mean
    ZO_af_subst{i}=ZO{i}-moyenne2{i};
    
    %Apply the SVD 
    [U{i}, S{i}, V{i}] = mySVD(ZO_af_subst{i});
    U_bar{i}=U{i}(:,1:40);   % Select the number of culomns that comprises the powred features (40). 
end
m=1;
i=1;
while i<=nb_representations
      if i<=nb_im_per_user
         ZO_Polaire_Reduced{i}=U_bar{m}'*ZO_Polaire{i};
         SI=size(ZO_Polaire_Reduced{i});
         Vector_Reduced{i} = reshape(ZO_Polaire_Reduced{i},1,SI(1)*SI(2));
         i=i+1;
      else
         m=m+1;
         nb_im_per_user=nb_im_per_user+8;
      end
end

end

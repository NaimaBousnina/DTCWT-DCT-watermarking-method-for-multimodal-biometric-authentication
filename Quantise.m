function Binary_Minutia_all = Quantise(Norm_Template)

%Number of vectors(=Number of users)
nb_vectors=size(Norm_Template);
for i=1:nb_vectors(2)
    nb_elements=length(Norm_Template{i}); %  le nombre des elements pour chaque vecteur.
end
for k=1:nb_vectors(2)   
    for i=1:nb_elements
        if Norm_Template{k}(i)>0
           Binary_Minutia_all{k}(i)=1;
        else
           Binary_Minutia_all{k}(i)=0;
        end
    end
end

end

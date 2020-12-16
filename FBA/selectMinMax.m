function NEWminmax = selectMinMax(Carbon,Nitrogen,Sulphur)

% mm1 = the first set of minmax (nitrogen)
% mm2 = the second set of minmax (carbon)
% bericht = vector, where each row is:
% Contains X: rxn with a met with and X formula
% balanced: same number of molecules left and right
% balanced - drain: exchange bewteen the cellualr and the extracellualr
% enviroment
% noflux: 0 flux
% set_bnd: LB and UB are defined
% unbalanced: different number of atoms
% unbalanced - CofPair: contains a cofactor
% no C/N in the reaction: this specific reaction doesn't contain any of the
% particualr atoms

mm1 = Nitrogen.newMinMax; 
mm2 = Carbon.newMinMax; 
mm3 = Sulphur.newMinMax; 
bericht1 = Nitrogen.report.summary;
bericht2 = Carbon.report.summary;
bericht3 = Sulphur.report.summary;
BoundsADJ1 = Nitrogen.bndAjusted.rawData;
BoundsADJ2 = Carbon.bndAjusted.rawData;
BoundsADJ3 = Sulphur.bndAjusted.rawData;
s1 = size(mm1,1);
s2 = size(mm1,1);
s3 = size(mm3,1);
if (s1 ~= s2 ) || (s1 ~= s3) || (s2~=s3)
    error('The minmax have different sizes; don''t know what that is possible')
end

for i=1:s1
    
    if ((BoundsADJ1(i,1) + BoundsADJ2(i,1) + BoundsADJ3(i,1) + BoundsADJ1(i,2) + BoundsADJ2(i,2) + BoundsADJ3(i,2)) == 0)
        NEWminmax(i,3) = {'FVA Bound'};
        NEWminmax(i,1) = {mm1(i,1)};
        NEWminmax(i,2) = {mm1(i,2)};
    end
   
    if ((isequal(bericht1{i,1}, 'Biomass nonproducing')) || (isequal(bericht1{i,1}, 'Biomass producing')) || (isequal(bericht1{i,1}, 'IgG Heavy Chain Synthesis')) || (isequal(bericht1{i,1}, 'IgG Light Chain Synthesis')))
        NEWminmax(i,1) = {mm1(i,1)};
        NEWminmax(i,2) = {mm1(i,2)};
        NEWminmax(i,3) = {'Important reaction'};

    elseif (isequal(bericht1{i,1}, 'Contains X') || isequal(bericht1{i,1}, 'balanced') || isequal(bericht1{i,1}, 'balanced - CofPair')|| isequal(bericht1{i,1}, 'balanced - drain'))
        % Set the stricter bounds
                
        if (mm1(i,1) < mm2(i,1)) && (mm1(i,1) < mm3(i,1))% Lower
            if mm2(i,1) < mm3(i,1)
                NEWminmax(i,1) = {mm3(i,1)};
                NEWminmax(i,3) = {'Sulphur'};
            else
                NEWminmax(i,1) = {mm2(i,1)};
                NEWminmax(i,3) = {'Carbon'};
            end
        elseif (mm2(i,1) < mm1(i,1)) && (mm2(i,1) < mm3(i,1))
            if mm1(i,1) < mm3(i,1)
                NEWminmax(i,1) = {mm3(i,1)};
                NEWminmax(i,3) = {'Sulphur'};
            else
                NEWminmax(i,1) = {mm1(i,1)};
                NEWminmax(i,3) = {'Nitrogen'};
            end
        else
            if mm1(i,1) < mm2(i,1)
                NEWminmax(i,1) = {mm2(i,1)};
                NEWminmax(i,3) = {'Carbon'};
            else
                NEWminmax(i,1) = {mm1(i,1)};
                NEWminmax(i,3) = {'Nitrogen'};
            end
        end
        if (mm1(i,2) > mm2(i,2)) && (mm1(i,2) > mm3(i,2)) % Upper
            if mm3(i,2) > mm2(i,2)
                NEWminmax(i,2) = {mm2(i,2)};
                NEWminmax(i,3) = {'Carbon'};
            else
                NEWminmax(i,2) = {mm3(i,2)};
                NEWminmax(i,3) = {'Sulphur'};
            end
        elseif (mm2(i,2) > mm1(i,2)) && (mm2(i,2) > mm3(i,2))
            if mm3(i,2) > mm1(i,2)
                NEWminmax(i,2) = {mm1(i,2)};
                NEWminmax(i,3) = {'Nitrogen'};
            else
                NEWminmax(i,2) = {mm3(i,2)};
                NEWminmax(i,3) = {'Sulphur'};
            end
        else
            if mm2(i,2) > mm1(i,2)
                NEWminmax(i,2) = {mm1(i,2)};
                NEWminmax(i,3) = {'Nitrogen'};
            else
                NEWminmax(i,2) = {mm2(i,2)};
                NEWminmax(i,3) = {'Carbon'};
            end
        end
    elseif (isequal(bericht1{i,1}, 'no N in the reaction') && isequal(bericht2{i,1}, 'no C in the reaction')) && (isequal(bericht3{i,1}, 'no S in the reaction'))
        NEWminmax(i,1) = {mm1(i,1)};
        NEWminmax(i,2) = {mm1(i,2)};
        NEWminmax(i,3) = {'No N and C and S'};
    elseif (isequal(bericht2{i,1}, 'no C in the reaction')) &&  (isequal(bericht1{i,1}, 'no N in the reaction'))
        NEWminmax(i,1) = {mm3(i,1)};
        NEWminmax(i,2) = {mm3(i,2)};
        NEWminmax(i,3) = {'No C and N'};  
    elseif (isequal(bericht2{i,1}, 'no C in the reaction')) &&  (isequal(bericht1{i,1}, 'no S in the reaction'))
        NEWminmax(i,1) = {mm1(i,1)};
        NEWminmax(i,2) = {mm1(i,2)};
        NEWminmax(i,3) = {'No C and S'};  
    elseif (isequal(bericht2{i,1}, 'no N in the reaction')) &&  (isequal(bericht1{i,1}, 'no S in the reaction'))
        NEWminmax(i,1) = {mm2(i,1)};
        NEWminmax(i,2) = {mm2(i,2)};
        NEWminmax(i,3) = {'No N and S'}; 
    elseif (isequal(bericht1{i,1}, 'no N in the reaction'))
        if mm2(i,1) < mm3(i,1)
            NEWminmax(i,1) = {mm3(i,1)};
        else
            NEWminmax(i,1) = {mm2(i,1)};
        end
        if mm2(i,2) < mm3(i,2)
            NEWminmax(i,2) = {mm2(i,2)};
        else
            NEWminmax(i,2) = {mm3(i,2)};
        end
        NEWminmax(i,3) = {'No N'};    
    elseif (isequal(bericht2{i,1}, 'no C in the reaction'))
        if mm1(i,1) < mm3(i,1)
            NEWminmax(i,1) = {mm3(i,1)};
        else
            NEWminmax(i,1) = {mm1(i,1)};
        end
        if mm1(i,2) < mm3(i,2)
            NEWminmax(i,2) = {mm1(i,2)};
        else
            NEWminmax(i,2) = {mm3(i,2)};
        end
        NEWminmax(i,3) = {'No C'};    
    elseif (isequal(bericht3{i,1}, 'no S in the reaction'))  
        if mm1(i,1) < mm2(i,1)
            NEWminmax(i,1) = {mm2(i,1)};
        else
            NEWminmax(i,1) = {mm1(i,1)};
        end
        if mm1(i,2) < mm2(i,2)
            NEWminmax(i,2) = {mm1(i,2)};
        else
            NEWminmax(i,2) = {mm2(i,2)};
        end
        NEWminmax(i,3) = {'No S'};    
    elseif (isequal(bericht1{i,1}, 'no_flux')|| isequal(bericht1{i,1}, 'set_bnd') ) 
        NEWminmax(i,1) = {mm1(i,1)};
        NEWminmax(i,2) = {mm1(i,2)};
        NEWminmax(i,3) = {'No bounds'}; 
    elseif (isequal(bericht1{i,1}, 'unbalanced') || isequal(bericht1{i,1}, 'unbalanced - CofPair')) 
        % Set the more loose bounds
        if (mm1(i,1) > mm2(i,1)) && (mm1(i,1) > mm3(i,1))% Lower
            if mm2(i,1) > mm3(i,1)
                NEWminmax(i,1) = {mm3(i,1)};
                NEWminmax(i,3) = {'Sulphur'};
            else
                NEWminmax(i,1) = {mm2(i,1)};
                NEWminmax(i,3) = {'Carbon'};
            end
        elseif (mm2(i,1) > mm1(i,1)) && (mm2(i,1) > mm3(i,1))
            if mm1(i,1) > mm3(i,1)
                NEWminmax(i,1) = {mm3(i,1)};
                NEWminmax(i,3) = {'Sulphur'};
            else
                NEWminmax(i,1) = {mm1(i,1)};
                NEWminmax(i,3) = {'Nitrogen'};
            end
        else
            if mm1(i,1) > mm2(i,1)
                NEWminmax(i,1) = {mm2(i,1)};
                NEWminmax(i,3) = {'Carbon'};
            else
                NEWminmax(i,1) = {mm1(i,1)};
                NEWminmax(i,3) = {'Nitrogen'};
            end
        end
        if (mm1(i,2) < mm2(i,2)) && (mm1(i,2) < mm3(i,2)) % Upper
            if mm3(i,2) < mm2(i,2)
                NEWminmax(i,2) = {mm2(i,2)};
                NEWminmax(i,3) = {'Carbon'};
            else
                NEWminmax(i,2) = {mm3(i,2)};
                NEWminmax(i,3) = {'Sulphur'};
            end
        elseif (mm2(i,2) < mm1(i,2)) && (mm2(i,2) < mm3(i,2))
            if mm3(i,2) < mm1(i,2)
                NEWminmax(i,2) = {mm1(i,2)};
                NEWminmax(i,3) = {'Nitrogen'};
            else
                NEWminmax(i,2) = {mm3(i,2)};
                NEWminmax(i,3) = {'Sulphur'};
            end
        else
            if mm2(i,2) < mm1(i,2)
                NEWminmax(i,2) = {mm1(i,2)};
                NEWminmax(i,3) = {'Nitrogen'};
            else
                NEWminmax(i,2) = {mm2(i,2)};
                NEWminmax(i,3) = {'Carbon'};
            end
        end
        NEWminmax(i,3) = {'Unbalanced'};   
    elseif (isequal(bericht1{i,1},'unconstr - oxPhos'))
        NEWminmax(i,1) = {-100};
        NEWminmax(i,2) = {100};
        NEWminmax{i,3} = {'Oxidative Phosporylation'}; 
    end
end

end

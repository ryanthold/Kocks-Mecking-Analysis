%% This is a function to convert Engineering Stress vs Strain to True Stress Strain
%% Inputs:
% ESS: x-by-2 double array double array of engineering stress[MPa]-strain[mm/mm] data with strain in the first column and stress in the second column
%% Outputs:
% TSS: x-by-2 double array double array of true stress[MPa]-strain[mm/mm] data with strain in the first column and stress in the second column
function[TSS]=ESS_to_TSS(ESS)
eng_strain=ESS(:,1);
eng_stress=ESS(:,2);

%Convert engineering stress and strain to true stress (ksi) and strain (%)
true_strain=log(eng_strain+1);
true_stress=eng_stress.*(1+eng_strain);

TSS=[true_strain,true_stress];
end
%% Code for the paper "Quantitative assessment of flow between cerebrospinal and interstitial fluid compartments in humans".
%
%  Code-authors: Anders Wåhlin, Viktor Vigren Näslund, Anders Eklund
%
%  gFlow is the function for estimating k1 and k2 from which CSF-to-ISF flow
%  rate and ISF volume fraction can be calculated.

%  getEstCt is a support function that calculates brain tissue concentrations given
%  a matrix of CSF concentrations (Nrois X time), as well as a given k1 and
%  k2 (or just as a function of k1 if extracellular volume fraction is fixed, using function suffix "_constantVe").
%
%  Implementation below is all intrathecal Gd flow calculations from the paper.

%% Main section
clear all
load dataset1.mat

for pat = 1:14; %loop cases
    for tissueClass = 1:3 %loop tissue classes


        if tissueClass ==1 %Cortex, estimate ve
            it = dataset1(pat).it;
            concCSF = dataset1(pat).concCsf_Ctx;
            concTissue = dataset1(pat).concCtx;
            dt = dataset1(pat).dt;
            points=dataset1(pat).points;
            tissueVolume = dataset1(pat).volCtx
            k0 = zeros(1,2)
            fun1=@(k)ITmodel(k,it,concCSF,concTissue,dt,points); %fit the model
            [kest,fval,exitflag,output] = fminunc(fun1,k0);
            q(pat,tissueClass) = kest(1)*tissueVolume %Glymphatic flow rate
            ve(pat,tissueClass) = kest(1)/kest(2) %ISF volume fraction

        end

        if tissueClass ==2 %White matter, keep ve constant
            it = dataset1(pat).it;
            concCSF = dataset1(pat).concCsf_Wm;
            concTissue = dataset1(pat).concWm;
            dt = dataset1(pat).dt;
            points=dataset1(pat).points;
            tissueVolume = dataset1(pat).volWm
            k0 = zeros(1)
            fun1=@(k)ITmodel(k,it,concCSF,concTissue,dt,points,ve(pat,1)); %fit the model
            [kest,fval,exitflag,output] = fminunc(fun1,k0);
            q(pat,tissueClass) = kest(1)*tissueVolume %Glymphatic flow rate
            ve(pat,tissueClass) = ve(pat,1) %ISF volume fraction
        end

        if tissueClass ==3 %Subcortical, keep ve constant
            it = dataset1(pat).it;
            concCSF = dataset1(pat).concCsf_Subcor;
            concTissue = dataset1(pat).concSubcor;
            dt = dataset1(pat).dt;
            points=dataset1(pat).points;
            tissueVolume = dataset1(pat).volSubcor
            k0 = zeros(1)
            fun1=@(k)ITmodel(k,it,concCSF,concTissue,dt,points,ve(pat,1)); %fit the model
            [kest,fval,exitflag,output] = fminunc(fun1,k0);
            q(pat,tissueClass) = kest(1)*tissueVolume %Glymphatic flow rate
            ve(pat,tissueClass) = ve(pat,1) %ISF volume fraction
        end

    end

    
end


%% Support functions

%%
function estCt = getEstCt(k,it,ccsf,dt)

% Estimates brain tissue concentrations provided a CSF concentrations
% (Nrois X time), k1 and k2.
% dt is the temporal resolution in hrs

for i = 1:size(ccsf,1) %loop over ROIs
    estCt(i,:) = k(1)*convolution(ccsf(i,:)',exp(-k(2).*it'))*dt ;
end

end
%%
function estCt = getEstCt_constantVe(k,it,ccsf,dt,ve)

% Estimates brain tissue concentrations provided a CSF concentrations
% (Nrois X time), k1 and k2.
% dt is the temporal resolution in hrs

for i = 1:size(ccsf,1) %loop over ROIs
    estCt(i,:) = k(1)*convolution(ccsf(i,:)',exp(-k(1)/ve.*it'))*dt ;
end

end

%%
function sse = ITmodel(k,it,ccsf,ct,dt,points,ve)
% Calculates the difference between measured (ct) and estimated brain tissue
% concentrations provided CSF concentrations (ccsf), k1 and k2 (stored as a
% vector k).
% The vector points specifies time points where the difference is calculated, and dt is the temporal resolution in hrs
switch nargin

    case 6 %ve estimated

        estCt = getEstCt(k,it,ccsf,dt);
        sse = sum(sum(((ct(:,points)-estCt(:,points)).^2)));

    case 7 %ve supplied

        estCt = getEstCt_constantVe(k,it,ccsf,dt,ve);
        sse = sum(sum(((ct(:,points)-estCt(:,points)).^2)));
end

end

%%
function c = convolution(a, b)
c = conv2(a(:), b(:), 'full');
c = c(1:size(a, 1));

if isrow(a)
    c = c.';
end
end

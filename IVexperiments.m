%% Code for the paper "Quantitative assessment of flow between cerebrospinal and interstitial fluid compartments in humans".
%
%  Code-authors: Anders Wåhlin, Viktor Vigren Näslund, Anders Eklund
%
%  gFlow is the function for estimating k1 and k2 from which CSF-to-ISF flow
%  rate and ISF volume fraction can be calculated.

%  getEstCt is a support function that calculates CSF  concentrations given
%  a matrix of tissue concentrations (Nrois X time), as well as a given k1 and
%  ve.
%
%  Implementation below describes intravenous Gd flow calculations from the paper.

%% Main section

load dataset2.mat 

for pat = 1:55; %loop cases


    it = dataset2(pat).it;
    concCSF = repmat(dataset2(pat).concCsf,2,1); %csf concentration. Since model expects a matrix, a dupliate row is added.
    concTissue = repmat(dataset2(pat).concTissue,2,1); %tissue concentration. Since model expects a matrix, a dupliate row is added.
    dt = dataset2(pat).dt;
    points = dataset2(pat).points;
    sasCsfVolume = dataset2(pat).vol

    k0 = zeros(1,2);
    fun1=@(k)IVmodel(k,it,concCSF,concTissue,dt,points); %fit the model
    [params,exitflag,output] = fmincon(fun1,k0,[],[],[],[],[-1 .15],[1 .35]);
    q(pat) = params(1)*sasCsfVolume
    ve(pat) = params(2);
    
    concTissueSum(pat) = sum(concTissue(1,:));
    concISFmeas = [0 concTissue(1,points)*(1/ve(pat))];
    concCSFmeas = [0 concCSF(1,points)];
    timeMeas = [0 it(points)]
    meanConcISF(pat,:) = mean(concTissue(1,:)*(1/ve(pat)));
    


end

idxAboveMedian = find(meanConcISF>median(meanConcISF));

figure(1)
histogram(q,-200:25:400)
figure(2)
histogram(q(idxAboveMedian),0:25:300)
mean(q(idxAboveMedian))

%%
function estCCSF = getEstCCSF(k,it,ct,dt)

% Estimates csf  concentrations provided tissue concentrations (Nrois X time), k1 and k2.
% dt is the temporal resolution in hrs

for i = 1:size(ct,1) %loop over ROIs
    estCCSF(i,:) = k(1)*convolution((ct(i,:))'*(1/k(2)),exp(-(k(1)).*it'))*dt;
end

end


%%
function sse = IVmodel(k,it,ccsf,ct,dt,points)

estCCSF = getEstCCSF(k,it,ct,dt);
sse = sum(sum((ccsf(:,points)-estCCSF(:,points)).^2))

end

%%
function c = convolution(a, b)
c = conv2(a(:), b(:), 'full');
c = c(1:size(a, 1));

if isrow(a)
    c = c.';
end
end

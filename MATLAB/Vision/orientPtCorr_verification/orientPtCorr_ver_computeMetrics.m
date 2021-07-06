function [corrSig,lineSig_1,lineSig_2] = orientPtCorr_ver_computeMetrics(x1,mid1,x2,mid2,F)

numElem = size(x1,1);
corrSig = zeros(numElem,1);
lineSig_1 = corrSig;
lineSig_2 = corrSig;

for k = 1:numElem
    corrSig(k,1) = abs([x2(k,:),1]*F*[x1(k,:)';1]);
    lineSig_1(k,1) = abs((x1(k,:)*[cos(mid1(k,2));sin(mid1(k,2))])-mid1(k,1));
    lineSig_2(k,1) = abs((x2(k,:)*[cos(mid2(k,2));sin(mid2(k,2))])-mid2(k,1));
end

%Whole DMN homogeneity

 function [WholeDMN,AnteriorDMN,posteriorDMN,APDMN]=DMN_Homogeneity(DMN,Ant_DMN,Post_DMN,Data)


%%%Imports
%DMN : DMN Mask
%Ant_DMN: Anterior DMN Mask
%Post_DMN: Posterior DMN Mask
%Data: 4-D volumetetric Data

% Exports
%r to z Correlations for
% WholeDMN: ALL of DMN
% AnteriorDMN: Anterior part of DMN
% PosteriorDMN: Posterior part of DMN
% APDMN: Anterior and Posterior part of DMN


%double type
Data=double(Data);

%Combine Anterior and Posterior Masks
AP_DMN=Ant_DMN+Post_DMN;


%Find Indices of Masks
 indsAll=find(DMN==1);
 nAll=length(indsAll);

 indsAnt=find(Ant_DMN==1);
 nAnt=length(indsAnt);

 indsPost=find(Post_DMN==1);
nPost=length(indsPost);

 indsAP=find(AP_DMN==1);
nAP=length(indsAP);

[ind1,ind2,ind3] = ind2sub(size(DMN),find(DMN == 1));
[ind1a,ind2a,ind3a] = ind2sub(size(Ant_DMN),find(Ant_DMN == 1));
[ind1p,ind2p,ind3p] = ind2sub(size(Post_DMN),find(Post_DMN == 1));
% [ind1ap,ind2ap,ind3ap] = ind2sub(size(AP_DMN),find(AP_DMN == 1));


 [s1,s2,s3,s4]=size(Data);


%Gets Time series at Indices Corresponding to Masks 
for i =1:nAll
 DataTrunc(i,:)=Data(ind1(i),ind2(i),ind3(i),:);
end

DataTruncPost=zeros(nPost,s4);
DataTruncAnt=zeros(nAnt,s4);

for ip =1:nPost
 DataTruncPost(ip,:)=Data(ind1p(ip),ind2p(ip),ind3p(ip),:);
end

for ia =1:nAnt
 DataTruncAnt(ia,:)=Data(ind1a(ia),ind2a(ia),ind3a(ia),:);
end

DataTruncAP=[DataTruncAnt;DataTruncPost];

%Computes Correlations and sums
% NHMat=zeros(size(DMN));
NHMatVals=(sum(corrcoef(DataTrunc'))-1)/(nAll-1);
NHMatValsAnt=(sum(corrcoef(DataTruncAnt'))-1)/(nAnt-1);
NHMatValsPost=(sum(corrcoef(DataTruncPost'))-1)/(nPost-1);
NHMatValsAP=(sum(corrcoef(DataTruncAP'))-1)/(nAP-1);  


%Fisher r to Z transformation and median
WholeDMN=median(atanh(NHMatVals));
AnteriorDMN=median(atanh(NHMatValsAnt));
posteriorDMN=median(atanh(NHMatValsPost));
APDMN=median(atanh(NHMatValsAP));

end

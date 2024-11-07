

function [AP_Assoc]=AP_Association(DMN,Ant_DMN,Post_DMN,Data)


%%%Imports
%DMN : DMN Mask
%Ant_DMN: Anterior DMN Mask
%Post_DMN: Posterior DMN Mask
%Data: 4-D volumetetric Data

% Exports
% AP_Assoc: Sum of correlations between voxels of Anterior and Posterior


%double type
Data=double(Data);

 indsAll=find(DMN==1);
 nAll=length(indsAll);

 indsAnt=find(Ant_DMN==1);
 nAnt=length(indsAnt);
 indsPost=find(Post_DMN==1);
nPost=length(indsPost);
NHValAllAnt_Post=[];
NHValAllPost_Ant=[];


 [ind1,ind2,ind3] = ind2sub(size(DMN),find(DMN == 1));
[ind1a,ind2a,ind3a] = ind2sub(size(Ant_DMN),find(Ant_DMN == 1));
[ind1p,ind2p,ind3p] = ind2sub(size(Post_DMN),find(Post_DMN == 1));

 [s1,s2,s3,s4]=size(Data);


DataTruncPost=zeros(nPost,s4);
DataTruncAnt=zeros(nAnt,s4);
tic





%Gets Posterior and Anterior Time Series
for ip =1:nPost
 DataTruncPost(ip,:)=Data(ind1p(ip),ind2p(ip),ind3p(ip),:);
end

for ia =1:nAnt
 DataTruncAnt(ia,:)=Data(ind1a(ia),ind2a(ia),ind3a(ia),:);
end


%Correlations cross-region between Anterior and Posterior
for ina =1:nAnt
    DataPostTemp=[DataTruncPost;DataTruncAnt(ina,:)];
    NHMatValsAnt_Post=(sum(corrcoef(DataPostTemp'))-1)/(nAnt);
    zCorrVal=atanh(NHMatValsAnt_Post(end));
    NHValAllAnt_Post(ina)= zCorrVal;
    {toc 100*ina/(nAnt+nPost) 'Ant'}
end

 for inp =1:nPost
    DataAntTemp=[DataTruncAnt;DataTruncPost(inp,:)];
    NHMatValsPost_Ant=(sum(corrcoef(DataAntTemp'))-1)/(nPost);
    zCorrVal=atanh(NHMatValsPost_Ant(end));
    NHValAllPost_Ant(inp)=zCorrVal;
    {toc 100*(inp+nAnt)/(nAnt+nPost) 'Post'}
 end

AP_Assoc=median([NHValAllPost_Ant,NHValAllAnt_Post]);


end


 



 

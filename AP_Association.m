

%NH
DMN=niftiread('C:\Users\leelab\Desktop\Conn_ICA\homo_output\DMN.nii'); %Whole DMN Mask
Ant_DMN=niftiread('C:\Users\leelab\Desktop\Conn_ICA\homo_output\anterior_DMN.nii'); %Anterior DMN Mask
Post_DMN=niftiread('C:\Users\leelab\Desktop\Conn_ICA\homo_output\posterior_DMN.nii'); % Posterior DMN Mask

 indsAll=find(DMN==1);
 nAll=length(indsAll);


 indsAnt=find(Ant_DMN==1);
 nAnt=length(indsAnt);

 indsPost=find(Post_DMN==1);
nPost=length(indsPost);

% totN=91*109*91;
% NHValAll=[];
NHValAllAnt_Post=[];
NHValAllPost_Ant=[];



 % [ind1,ind2,ind3] = ind2sub(size(DMN),find(DMN == 1));
[ind1a,ind2a,ind3a] = ind2sub(size(Ant_DMN),find(Ant_DMN == 1));
[ind1p,ind2p,ind3p] = ind2sub(size(Post_DMN),find(Post_DMN == 1));

tic
files=dir('*pre*.nii');
num=1;
for subj=1:length(files)
    subject=files(subj).name(1:10);
    


    for p=1:2
        
        NHMat=zeros(size(DMN));

        if p==1
            {subject subj 'pre' toc}
            subjfile=dir(strcat(subject,'*pre*.nii'));
        else
           {subject subj 'post' toc}

            subjfile=dir(strcat(subject,'*post*.nii'));
        end

 Data=double(niftiread(subjfile.name));
DataTruncPost=zeros(nPost,300);
DataTruncAnt=zeros(nAnt,300);

for ip =1:nPost
 DataTruncPost(ip,:)=Data(ind1p(ip),ind2p(ip),ind3p(ip),:);
end

for ia =1:nAnt
 DataTruncAnt(ia,:)=Data(ind1a(ia),ind2a(ia),ind3a(ia),:);
end

for ina =1:nAnt
    DataPostTemp=[DataTruncPost;DataTruncAnt(ina,:)];
    NHMatValsAnt_Post=(sum(corrcoef(DataPostTemp'))-1)/(nAnt);
    zCorrVal=atanh(NHMatValsAnt_Post(end));
    NHValAllAnt_Post(ina,num)= zCorrVal;
    NHMat(ind1a(ina),ind2a(ina),ind3a(ina))=zCorrVal;
     {subj p toc 100*ina/(nAnt+nPost) 'Ant'}
end


 for inp =1:nPost
    DataAntTemp=[DataTruncAnt;DataTruncPost(inp,:)];
    NHMatValsPost_Ant=(sum(corrcoef(DataAntTemp'))-1)/(nPost);
    zCorrVal=atanh(NHMatValsPost_Ant(end));
    NHValAllPost_Ant(inp,num)=zCorrVal;
    NHMat(ind1p(inp),ind2p(inp),ind3p(inp))=zCorrVal;
    {subj p toc 100*(inp+nAnt)/(nAnt+nPost) 'Post'}
 end


        if p==1
            outputStruct5.(subject).pre=NHMat;
        else
            outputStruct5.(subject).post=NHMat;
        end


num=num+1;

    end
end

AP_Assoc_All=median([NHValAllPost_Ant;NHValAllAnt_Post]);

AP_Assoc_Pre=APDMN(1:2:35)';
AP_Assoc_Post=APDMN(2:2:36)';

  save('outputStructTest5.mat','outputStruct5');



 

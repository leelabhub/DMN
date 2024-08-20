%Whole DMN homogeneity



DMN=niftiread('C:\Users\leelab\Desktop\Conn_ICA\homo_output\DMN.nii');
Ant_DMN=niftiread('C:\Users\leelab\Desktop\Conn_ICA\homo_output\anterior_DMN.nii');
Post_DMN=niftiread('C:\Users\leelab\Desktop\Conn_ICA\homo_output\posterior_DMN.nii');
AP_DMN=Ant_DMN+Post_DMN;


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
[ind1ap,ind2ap,ind3ap] = ind2sub(size(AP_DMN),find(AP_DMN == 1));



tic
files=dir('*pre*.nii');
num=1;
for subj=1:length(files)
    subject=files(subj).name(1:10);
    


    for p=1:2
        
        % NHMat=zeros(size(DMN));

        if p==1
            {subject subj 'pre' toc}
            subjfile=dir(strcat(subject,'*pre*.nii'));
        else
           {subject subj 'post' toc}

            subjfile=dir(strcat(subject,'*post*.nii'));
        end

 Data=double(niftiread(subjfile.name));
for i =1:nAll
 DataTrunc(i,:)=Data(ind1(i),ind2(i),ind3(i),:);
end

DataTruncPost=zeros(nPost,300);
DataTruncAnt=zeros(nAnt,300);

for ip =1:nPost
 DataTruncPost(ip,:)=Data(ind1p(ip),ind2p(ip),ind3p(ip),:);
end

for ia =1:nAnt
 DataTruncAnt(ia,:)=Data(ind1a(ia),ind2a(ia),ind3a(ia),:);
end

DataTruncAP=[DataTruncAnt;DataTruncPost];


% NHMat=zeros(size(DMN));
NHMatVals=(sum(corrcoef(DataTrunc'))-1)/(nAll-1);
NHValAll(:,num)=NHMatVals;

NHMatValsAnt=(sum(corrcoef(DataTruncAnt'))-1)/(nAnt-1);
NHValAllAnt(:,num)=NHMatValsAnt;


NHMatValsPost=(sum(corrcoef(DataTruncPost'))-1)/(nPost-1);
NHValAllPost(:,num)=NHMatValsPost;


NHMatValsAP=(sum(corrcoef(DataTruncAP'))-1)/(nAP-1);
NHValAllAP(:,num)=NHMatValsAP;


for i =1:nAnt
    NHMat(ind1a(i),ind2a(i),ind3a(i))=atanh(NHMatValsAnt(i));
end

for i =1:nPost
    NHMat(ind1p(i),ind2p(i),ind3p(i))=atanh(NHMatValsPost(i));
end

% % for i =1:nAP
% %     NHMat(ind1AP(i),ind2AP(i),ind3AP(i))=atanh(NHMatValsAP(i));
% % end



        if p==1
            outputStruct.(subject).pre=NHMat;
        else
            outputStruct.(subject).post=NHMat;
        end


num=num+1;

end
end


WholeDMN=mean(atanh(NHValAll));
WholeDMN_Pre=WholeDMN(1:2:35)';
WholeDMN_Post=WholeDMN(2:2:36)';



AnteriorDMN=mean(atanh(NHValAllAnt));
AnteriorDMN_Pre=AnteriorDMN(1:2:35)';
AnteriorDMN_Post=AnteriorDMN(2:2:36)';

posteriorDMN=mean(atanh(NHValAllPost));
posteriorDMN_Pre=posteriorDMN(1:2:35)';
posteriorDMN_Post=posteriorDMN(2:2:36)';


APDMN=mean(atanh(NHValAllAP));
APDMN_Pre=APDMN(1:2:35)';
APDMN_Post=APDMN(2:2:36)';


AnteriorDMN_med=median(atanh(NHValAllAnt));
AnteriorDMN_med_Pre=AnteriorDMN_med(1:2:35)';
AnteriorDMN_med_Post=AnteriorDMN_med(2:2:36)';

posteriorDMN_med=median(atanh(NHValAllPost));
posteriorDMN_med_Pre=posteriorDMN_med(1:2:35)';
posteriorDMN_med_Post=posteriorDMN_med(2:2:36)';


Ratio_Pre=AnteriorDMN_med_Pre./posteriorDMN_med_Pre);
Ratio_Post=AnteriorDMN_med_Post./posteriorDMN_med_Post;





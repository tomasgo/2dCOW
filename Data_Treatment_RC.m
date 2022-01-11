%%%%%%%%%%%%%%%%%%%%%%%%%%%  DATA TREATMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   REAL CASES SPILLS %%%%%%%%%%%%%%%%%%%%%%%%%

%% Info about the file

% Data.XalunN1 = entire row-wise normalized (norm = 1) data
% Data.XalunDSN1 = Data.XalunN1 after downscaling (* 1/RSDr)
% Data.ind = index delimitating each SIC
% Data.Xalun2N1 = Data containing only signficant SICs and row-wise
% normalized

% LocalPCA = local PCA models of each SIC to compare total variance with
% replicate variances


% WPCA = PCA on Data.XalunDSN1 (samples replicates taken as validation set)
% WPCA2 = PCA on Data.Xalun2N1
% WPCA2_eva = PCA considering evaporation




%% Selecting M/Z ions

ions = [113.15;127.15;183.2;197.25;83.1;81.05;95.10;109.1;137.15;123.10;163.15;91.05;...
    105.05;119.10;134.10;147.10;161.15;117.05;131.10;145.10;159.10;174.15;188.15;...
    129.05;143.10;157.10;171.10;185.15;167.10;182.10;128.05;142.10;156.10;170.10;...
    184.10;184.15;166.10;180.10;194.10;208.10;178.10;192.10;206.10;220.10;234.15;184.05;...
    198.05;212.05;226.10;202.10;216.10;230.10;228.10;74.05];

Samples.Batch1 = {'QC2.cdf';'QC1.cdf';'Blank1.cdf';'3211-4a.cdf';'3210-1.cdf';'3203-4.cdf';...
    '3203-3.cdf';'3188-2.cdf';'3184-2.cdf';'3184-1b.cdf';'3184-1a.cdf';'3183-2.cdf';...
    '3183-1.cdf';'3182-3.cdf';'3182-2b.cdf';'3182-2a.cdf';'3178-4.cdf'};

Samples.Batch2 = {'QC4.cdf';'QC3.cdf';'Blank2.cdf';'3243-8.cdf';'3243-7.cdf';'3243-4b.cdf';...
    '3243-4a.cdf';'3242-5.cdf';'3242-3b.cdf';'3242-3a.cdf';'3242-2.cdf';'3229-1.cdf';...
    '3221-2b.cdf';'3221-2a.cdf';'3221-1.cdf';'3211-5.cdf';'3211-4b.cdf'};

Samples.Batch3 = {'QC5.cdf';'J3261.cdf';'J3253-23b.cdf';'J3253-23a.cdf';'J3253-20b.cdf';...
    'J3253-20a.cdf';'J3253-19.cdf';'J3253-13.cdf';'J3253-12b.cdf';'J3253-12a.cdf';...
    'J3253-8.cdf';'Blank3.cdf';'11781.cdf';'3243-12b.cdf';'3243-12a.cdf';'3243-11.cdf'};


%% Loading the High Res GC-QTOF CDF data

for i = 1:length(Samples.Batch1);
    
    [Batch1.acqT{i}, Batch1.SIC{i}] = readHRGCSIC(Samples.Batch1{i},ions);
       
end

clear i

%% Combining SICs

sic_groups = {1:4;5;6:9;10;11;12:17;18:23;24:28;29:30;31;32;33;34;35:36;37;38;39;40;41;...
    42;43;44;45;46;47;48;49;50;51;52;53;54};

for i = 1:length(Samples.Batch1);
    Batch1.SIC_G{i} = zeros(size(Batch1.SIC{i},1),length(sic_groups));
    
    for j = 1:length(sic_groups);
        Batch1.SIC_G{i}(:,j) = sum(Batch1.SIC{i}(:,sic_groups{j}),2);   
    end
    
end
clear j i


%% Refolding TIC (summed SIC) to GCxGC for future phase-shifting
for i = 1:length(Samples.Batch3);
    
    [Batch1.TIC_2D{i},Batch1.rt1{i},Batch1.rt2{i},~] = reshape2DChrom(sum(Batch1.SIC_G{i},2),Batch1.acqT{i},7,[]);

end
clear i


%% Phase-Shifting TIC

i = 16;
Batch3.shift{i} = size(Batch3.TIC_2D{i},1)-89;
X = Batch3.TIC_2D{i};

Xn = circshift(X,Batch3.shift{i});
Xn2 = circshift(Xn(Batch3.shift{i}+1:end,:),-1,2);
Xn(Batch3.shift{i}+1:end,:) = Xn2;

figure;imagesc(Xn)
axis xy
clear X Xn2 i Xn
set(gcf,'Color',[1,1,1])

%% Phase-Shifting SICs

for i = 1:length(Samples.Batch3);

    for j = 1:length(sic_groups);
   
        [X_2D,~,~,~] = reshape2DChrom([Batch3.SIC_G{i}(length(Batch3.rt2{i})-Batch3.shift{i}:end,j);...
            Batch3.SIC_G{i}(size(Batch3.SIC_G{i},2)*ones(length(Batch3.rt2{i})-Batch3.shift{i}+1,1),j)]...
            ,Batch3.acqT{i},7,[]); 
        
        Batch3.SIC_GPS{i}(:,:,j) = X_2D;
        clear X_2D
    end
    
end

clear i j 

%% Checking for rt1 shifts
i = 1;
X(:,:,1) = sum(Xal{1,i},3);
X(:,:,2) = sum(Xal{2,i},3);
X(:,:,3) = sum(Xal{17,i},3);
X(:,:,4) = sum(Xal{32,i},3);
X(:,:,5) = Batch2.SIC_GPSsub{2,i};

for i = 1:size(X,3);
    X2(:,i) = sum(X(:,:,i),1,'omitnan');
end

figure;plot(normc(X2(:,1:5)))
axis tight
clear X X2 i


%% Setting the 2D-RT limits for 2D alignment

rtLim{1} = [1 350;1 167];
rtLim{2} = [1 350;1 167];
rtLim{3} = [1 320;20 167];
rtLim{4} = [50 320;30 140];
rtLim{5} = [85 270;10 80];
rtLim{6} = [1 325;1 120];
rtLim{7} = [70 300;1 80];
rtLim{8} = [80 320;1 80];
rtLim{9} = [200 340;1 70];
rtLim{10} = [80 320;1 60];
rtLim{11} = [175 195;1 30];
rtLim{12} = [195 250;1 40];
rtLim{13} = [210 280;1 60];
rtLim{14} = [225 315;10 60];
rtLim{15} = [250 260;10 30];
rtLim{16} = [215 360;10 60];
rtLim{17} = [240 380;10 70];
rtLim{18} = [250 380;15 70];
rtLim{19} = [270 310;15 50];
rtLim{20} = [310 335;20 50];
rtLim{21} = [325 360;30 70];
rtLim{22} = [330 380;30 80];
rtLim{23} = [350 400;35 80];
rtLim{24} = [290 305;15 40];
rtLim{25} = [305 330;20 50];
rtLim{26} = [320 350;25 60];
rtLim{27} = [330 365;30 70];
rtLim{28} = [360 370;30 50];
rtLim{29} = [365 395;35 75];
rtLim{30} = [380 420;40 90];
rtLim{31} = [415 425;50 75];
rtLim{32} = [230 320;50 140];


%% Selecting the sub-regions of the SIC according to rtLim for data handling

for i = 1:length(rtLim);
    
    for j = 1: length(Batch3.SIC_GPS);
    Batch3.SIC_GPSsub{j,i} = Batch3.SIC_GPS{j}(rtLim{i}(2,1):rtLim{i}(2,2),...
        rtLim{i}(1,1):rtLim{i}(1,2),i);
    end
    
end
clear i j

%% Replacing NaNs due to reshaping the 2D Chromatograms
A = Batch3.SIC_GPSsub;

for i = 1:3;
    
    for j = 1:size(A,1);
        [c] = find(isnan(A{j,i}(end,:)));
        
        A{j,i}(end,c) = A{j,i}(end-1,c);
    
    end

end

Batch3.SIC_GPSsub = A; clear A c i j



%% PEAK ALIGNMENT using 2D-COW for TESTING - using QC samples only

% Concatenating the images of each SIC from all the Batches
X = [];
X(:,:,1) = sum(Batch1.SIC_GPS{2},3);
X(:,:,2) = sum(Batch1.SIC_GPS{1},3);
X(:,:,3) = sum(Batch2.SIC_GPS{2},3);
X(:,:,4) = sum(Batch2.SIC_GPS{1},3);
X(:,:,5) = sum(Batch3.SIC_GPS{1},3);

X = X(1:end-1,:,:);
for i = 1:5;
[retX(:,:,i),~] = TwoDCOW(X(:,:,i),X(:,:,3),[50,5],[10,1]);


end

%% %%%%%%%%%%% PEAK ALIGNMENT of each 2D-SIC using 2D-COW %%%%%%%%%%%%%%%%

%% Defining the Segments and Slack Limits

sic_group = 32;

for i = 1:16;
    plot(reshape(Batch1.SIC_GPSsub{i,sic_group},size(Batch1.SIC_GPSsub{i,sic_group},1)*size(Batch1.SIC_GPSsub{i,sic_group},2),1))
    hold on
    plot(reshape(sum(Batch2.SIC_GPSsub{i,sic_group},3),size(Batch2.SIC_GPSsub{i,sic_group},1)*size(Batch2.SIC_GPSsub{i,sic_group},2),1))
    hold on
    plot(reshape(sum(Batch3.SIC_GPSsub{i,sic_group},3),size(Batch3.SIC_GPSsub{i,sic_group},1)*size(Batch3.SIC_GPSsub{i,sic_group},2),1))
    plot(reshape(sum(Batch2.SIC_GPSsub{2,sic_group},3),size(Batch2.SIC_GPSsub{2,sic_group},1)*size(Batch2.SIC_GPSsub{2,sic_group},2),1)...
        ,'--','Color','Black','LineWidth',2)
end

clear i
% Defining the limits for the warping optimization

COW2D.seglimy{sic_group} = 5;   % max segment length for y (rt1)
COW2D.slacky{sic_group} = 1;    % max slack for y (rt1)
COW2D.seglimx{sic_group} = 10:2:20; % 5 pts spaced segment length for x (rt2)
COW2D.slackx{sic_group} = 1:1:5;   % 1 pt spaced slack for x (rt1)

clear Xrefunf peakloc sic_group peakw

%% RUNNING 2D-COW with optimization of the paramters

sic_group = 32;

for m = 1:size(Batch1.SIC_GPSsub,1);    
    X1(:,:,m) = Batch1.SIC_GPSsub{m,sic_group};
end
X1n = X1(:,:,[1,2,4:size(X1,3)]);
X1 = X1n; clear X1n

for m = 1:size(Batch2.SIC_GPSsub,1);
    X2(:,:,m)= Batch2.SIC_GPSsub{m,sic_group};
end
   X2n = X2(:,:,[1,4:size(X2,3)]);
   X2 = X2n; clear X2n
   
for m = 1:size(Batch3.SIC_GPSsub,1);    
    X3(:,:,m) = Batch3.SIC_GPSsub{m,sic_group};
end
X3n = X3(:,:,[1:11,13:size(X3,3)]);
X3 = X3n; clear X3n

X = cat(3,X1,X2,X3);clear X1 X2 X3

%Activate only when running first SIC
%COW2D.WEff = cell(size(X,3),length(sic_groups));
%COW2D.S = cell(size(X,3),length(sic_groups));
%COW2D.Pf = cell(size(X,3),length(sic_groups));
%Xal = cell(length(Samples.Batch1)+length(Samples.Batch2)+length(Samples.Batch3),length(sic_groups));
% ----------------------------------------------------------------

Xref = Batch2.SIC_GPSsub{2,sic_group}; % QC3 is set as ref chrom
Xrefunf = reshape(Xref,1,size(Xref,1)*size(Xref,2));
Xwunf = zeros(2,size(Xrefunf,1)*size(Xrefunf,2));
Xwunf(2,:) = Xrefunf;

for i = 1:size(X,3);
normX = norm(reshape(X(:,:,i),1,size(X,1)*size(X,2)));
Xw = cell(length(COW2D.seglimx{sic_group}),length(COW2D.slackx{sic_group}));

    for j = 1:length(COW2D.seglimx{sic_group});
   
        for k = 1:length(COW2D.slackx{sic_group});
       
            inSLen = [COW2D.seglimx{sic_group}(j),COW2D.seglimy{sic_group}];
            inMaxW = [COW2D.slackx{sic_group}(k),COW2D.slacky{sic_group}];
            
            if COW2D.seglimx{sic_group}(j) <= COW2D.slackx{sic_group}(k);
               COW2D.S{i,sic_group}(j,k) = 0;
               COW2D.Pf{i,sic_group}(j,k) = 0;
            else
                [Xw{j,k},~] = TwoDCOW(X(:,:,i),Xref,inSLen,inMaxW);
                Xwunf(1,:) = reshape(Xw{j,k},1,size(Xw{j,k},1)*size(Xw{j,k},2));
            
                % Computing Simplicity
                COW2D.S{i,sic_group}(j,k) = (sum(svds(Xwunf/sqrt(sum(Xwunf(:).^2)),1).^4));
            
            % Computing Peak Factor
                COW2D.Pf{i,sic_group}(j,k) = abs(((norm(Xwunf(1,:))-normX))/normX);
                COW2D.Pf{i,sic_group}(j,k) = (1-min([COW2D.Pf{i,sic_group}(j,k) 1]))^2;
            
                clear inSLen inMaxW Xwunf
            end
        end
    
    end
    % Computing the Warping Effect
    COW2D.WEff{i,sic_group} = COW2D.S{i,sic_group} + COW2D.Pf{i,sic_group};
    
    mweff = max(COW2D.WEff{i,sic_group}(:));
    [r,c] = find(COW2D.WEff{i,sic_group} == mweff);
    
    Xal{i,sic_group} = Xw{r,c}; 
    
    clear r c mweff normX Xw
end
clear sic_group X Xref Xrefunf i j k m X

%% Creating the huge augmented matrix - NORMALIZING ENTIRE ROW-VECTOR

Data.Xalun = [];
for j = 1:32;
    for i = 1:46;
        X1(i,:) = reshape(Xal{i,j},1,size(Xal{i,j},1)*size(Xal{i,j},2));
    end
    X1(47,:) = reshape(Batch2.SIC_GPSsub{2,j},1,size(Batch2.SIC_GPSsub{2,j},1)*size(Batch2.SIC_GPSsub{2,j},2));
    Data.Xalun = [Data.Xalun,X1];
    clear X1
end
%Data.Xalun = normr(Data.Xalun);
clear i j

%% Inspecting SICs
k = 11;
figure;plot(Data.Xalun(:,Data.ind(k)+1:Data.ind(k+1))')

%%
ind(1:10) = ones(1,10).*29220;



%% REFINING PEAKS ALIGNMENTS USING 1D-ICOSHIFT

% FINDING PEAKS AND CREATING PEAK LIST

k = 27;
[peakinfo{k}, peakloc{k}] = findpeaks(Data.Xalun(47,Data.ind(k)+1:Data.ind(k+1)),'minpeakheight',5E4, ... 
    'minpeakdistance',10);

figure;plot(Data.Xalun(:,Data.ind(k)+1:Data.ind(k+1))')
hold on
plot(peakloc{k},peakinfo{k},'*')
axis tight



%% Creating intervals vector 
k = 32;

if length(peakinfo{k}) == 1;
    
    int{k}(1,1:2) = [1 length(Data.Xalun(47,Data.ind(k)+1:Data.ind(k+1)))];
    
else
    
    int{k}(1,1:2) = [1 round(mean([peakloc{k}(1,1);peakloc{k}(1,2)]))];
 
    for i = 2:length(peakloc{k})-1;
        int{k}(1,length(int{k})+1) = int{k}(1,length(int{k}))+1;
        int{k}(1,length(int{k})+1) = round(mean([peakloc{k}(1,i);peakloc{k}(1,(i)+1)]));
    end
 
    int{k}(1,length(int{k})+1) = int{k}(1,length(int{k}))+1;
    int{k}(1,length(int{k})+1) = Data.ind(1,(k+1)) - Data.ind(1,k);  % add +1 for k = 1
end

%% Running 1D-COW in each interval

k = 32;

X1 = Data.Xalun(1:46,Data.ind(k)+1:Data.ind(k+1));
Xref = Data.Xalun(47,Data.ind(k)+1:Data.ind(k+1));

[~,~,al2] = optim_cow(X1(:,7000:7600),[20,50,5,10],[0 3 50 0.15],Xref(:,7000:7600));

clear X1 Xref k

%% RUNNING iCOSHIFT

k = 32;

[xCS{k},~,~] = icoshift(Data.Xalun(47,Data.ind(k)+1:Data.ind(k+1)),Data.Xalun(1:46,Data.ind(k)+1:Data.ind(k+1)),int{k},'b');


figure;
subplot(2,1,2); plot(xCS{k}'); axis tight
subplot(2,1,1); plot(Data.Xalun(:,Data.ind(k)+1:Data.ind(k+1))'); axis tight



%% Concatenating the RE-ALIGNED chromatograms into a single Augmented Matrix

Data.Xalun2 = [];

for i = 1:length(xCS)
   
    Data.Xalun2 = [Data.Xalun2 xCS{i}];
    
end

Data.Xalun2(47,:) = Data.Xalun(47,:);


%% REFOLDING THE RE-ALIGNED CHROMATOGRAMS TO THE ORIGINAL 2D-structure

for i = 1:length(WPCA.sics);
    
    for j = 1:size(Data.Xalun2,1);
       
        Xal2{j,i} = reshape(Data.Xalun2(j,Data.ind(WPCA.sics(i))+1:Data.ind(WPCA.sics(i)+1)),...
            size(Xal{1,WPCA.sics(i)},1),size(Xal{1,WPCA.sics(i)},2));
        
    end
    
    
end



%% Creating the huge augmented matrix - Normalizing each SIC individually

Data.XalunN2 = [];
for j = 1:32;
    for i = 1:46;
        X1(i,:) = reshape(Xal{i,j},1,size(Xal{i,j},1)*size(Xal{i,j},2));
    end
    X1(47,:) = reshape(Batch2.SIC_GPSsub{2,j},1,size(Batch2.SIC_GPSsub{2,j},1)*size(Batch2.SIC_GPSsub{2,j},2));
    X1 = normr(X1);
    Data.XalunN2 = [Data.XalunN2,X1];
    clear X1
end
clear i j

%% Inspecting aligned QC samples

QC = ([Data.Xalunb(2,Data.ind(16):Data.ind(17));Data.Xalunb(1,Data.ind(16):Data.ind(17));Data.Xalunb(17,Data.ind(16):Data.ind(17));Data.Xalunb(47,Data.ind(16):Data.ind(17));Data.Xalunb(32,Data.ind(16):Data.ind(17))]);

figure;plot(QC')

clear QC


%% Creating Dataset

s = [Samples.Batch1([1,2,4:end],1);Samples.Batch2([1,4:end],1);Samples.Batch3([1:11,13:end],1);Samples.Batch2(2,1)];
%%
%editds(B);
Data.label{1,1} = strtok(s,'.');
Data.class{1,1} = class;

%% Extracting the peaks for WLS-PCA

% Establishing the threshold for the peaks

figure;findpeaks(mean(Data.Xalun),'MinPeakHeight',1E-3,'MinPeakWidth',2,'Annotate','extent')

%% Extracting the peaks

[peaks.int,peaks.loc,peaks.w,~] = findpeaks(mean(Data.Xalun),'MinPeakHeight',1E-3,'MinPeakWidth',2);


%% Extracting only variables belonging to the main peaks

peaks.w2 = zeros(1,length(peaks.w));
for i = 2:length(peaks.loc)-1;
       
    if peaks.loc(i)-floor(peaks.w(i)/2)<= peaks.loc(i-1)+floor(peaks.w(i-1)/2);
       X1 = diff(mean(Data.Xalun(:,floor(mean([peaks.loc(i-1),peaks.loc(i)])):floor(mean([peaks.loc(i),peaks.loc(i)+peaks.w(i)/2])))));
       [~,a1] = find(X1 >= 0,1,'first');
       [~,a2] = find(X1 <= 0,1,'last');
       if isempty(a2) == 1;
           [~,a2] = find(X1 >= 0,1,'last');
       end
       peaks.w2(i) = a2-a1;
       clear X1 a1 a2
     
    elseif peaks.loc(i)+floor(peaks.w(i)/2)>= peaks.loc(i+1)-floor(peaks.w(i+1)/2);
        X1 = diff(mean(Data.Xalun(:,floor(mean([peaks.loc(i)-peaks.w(i)/2,peaks.loc(i)])):floor(mean([peaks.loc(i),peaks.loc(i+1)])))));
        [~,a1] = find(X1 >= 0,1,'first');
        [~,a2] = find(X1 >= 0,1,'last');
        peaks.w2(i) = a2-a1;
       clear X1 a1 a2  
       
     elseif peaks.loc(i)-floor(peaks.w(i)/2)<= peaks.loc(i-1)+floor(peaks.w(i-1)/2) && peaks.loc(i)+floor(peaks.w(i)/2)>= peaks.loc(i+1)-floor(peaks.w(i+1)/2);
       X1 = diff(mean(Data.Xalun(:,floor(mean([peaks.loc(i-1),peaks.loc(i)])):floor(mean([peaks.loc(i),peaks.loc(i+1)])))));
       [~,a1] = find(X1 >= 0,1,'first');
       [~,a2] = find(X1 >= 0,1,'last');
       peaks.w2(i) = a2-a1;
       clear X1 a1 a2
    end
    
end

[~,a1] = find(peaks.w2);
peaks.w(a1) = peaks.w2(a1);

clear a1

figure;plot(mean(Data.Xalun)),hold on,herrorbar(peaks.loc,peaks.int,peaks.w/2,'*')


%% Extracting the signals

Data.XalunB = [];

for i = 1:length(peaks.loc);
   
    X1 = Data.Xalun(:,peaks.loc(i)-floor(peaks.w(i)/2):peaks.loc(i)+floor(peaks.w(i)/2));
    
    Data.XalunB = [Data.XalunB X1];
    
end

%% Calculating the Weights for WLS-PCA

X1 = ([Data.Xalun2(2,:);Data.Xalun2(1,:);Data.Xalun2(17,:);Data.Xalun2(47,:);Data.Xalun2(32,:)]);
X2 = ([Data.Xalun2(3,:);Data.Xalun2(31,:)]);
X3 = ([Data.Xalun2(9,:);Data.Xalun2(10,:)]);
X4 = ([Data.Xalun2(14,:);Data.Xalun2(15,:)]);
X5 = ([Data.Xalun2(20,:);Data.Xalun2(21,:)]);
X6 = ([Data.Xalun2(23,:);Data.Xalun2(24,:)]);
X7 = ([Data.Xalun2(27,:);Data.Xalun2(28,:)]);
X8 = ([Data.Xalun2(34,:);Data.Xalun2(35,:)]);
X9 = ([Data.Xalun2(36,:);Data.Xalun2(37,:)]);
X10 = ([Data.Xalun2(40,:);Data.Xalun2(41,:)]);
X11 = ([Data.Xalun2(44,:);Data.Xalun2(45,:)]);


A = sqrt((((size(X1,1)-1)*var(X1))+((size(X2,1)-1)*var(X2))+((size(X3,1)-1)*var(X3))+...
    ((size(X4,1)-1)*var(X4))+((size(X5,1)-1)*var(X5))+((size(X6,1)-1)*var(X6))+((size(X7,1)-1)*var(X7))...
    +((size(X8,1)-1)*var(X8))+((size(X9,1)-1)*var(X9))+((size(X10,1)-1)*var(X10))+((size(X11,1)-1)*var(X11)))./...
    (size(X1,1)+size(X2,1)+size(X3,1)+size(X4,1)+size(X5,1)+size(X6,1)+size(X7,1)+size(X8,1)+size(X9,1)+size(X10,1)+size(X11,1)-11));

clear X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11

%W = std(Xalun2)./A;figure;plot(W')
WPCA2B.W = 1./A;figure;plot(WPCA2B.W')
clear A



%% Calculating the Weights for WLS-PCA (ONLY CONSIDERING QCs)

X1 = ([Data.Xalun2N1(2,:);Data.Xalun2N1(1,:);Data.Xalun2N1(17,:);Data.Xalun2N1(47,:);Data.Xalun2N1(32,:)]);

LocalPCA.W = 1./(std(X1)./mean(X1));figure;plot(LocalPCA.W')

clear X1

%% Downscaling the variables and replacing NaNs when necessary

X = zeros(size(Data.Xalun2,1),size(Data.Xalun2,2));
for i = 1:size(Data.Xalun2,1);
    X(i,:) = Data.Xalun2(i,:).*LocalPCA.W;
end

a = find(isnan(LocalPCA.W));
if isempty(a) == 0;
    for i = 1:size(X,1);
        X(i,a) = zeros(1,length(a));
    end
end
clear a

%plot(X')

Data.Xalun2DS = X;  clear X

%% Defining Calibration and Validation Sets - PCA after simple downscaling 

editds(Data.Xalun2DS);
%%
Data.Xalun2DS = DataXalun2DS;  clear DataXalun2DS,

Data.Xalun2DS.label{1,1} = strtok(Data.info(:,1),'.');
Data.Xalun2DS.class{1,1} = Data.info(:,2);


%% Defining Calibration and Validation Sets - PCA after simple downscaling

LocalPCA.val = find(cell2mat(Data.info(:,3))); %Samples in Data.info(:,3) extracted for validation

LocalPCA.cal = find(cell2mat(Data.info(:,3)) == 0);

%%  Running WPCA (GUI)

Xcal = Data.XalunDS(WPCA2B.cal,:);
Xval = Data.XalunDS(WPCA2B.val,:);

pca

%% Plotting the WPC1 x WPC2 scores (with validation)


qc = cellfun('isempty', strfind(Data.info(WPCA2.cal,2),'QC'));
qc2 = find(qc == 0);
A = Data.info(WPCA2.cal,1);
figure;plot(WPCA2model.loads{1,1}(qc2,1),WPCA2model.loads{1,1}(qc2,2),'o','MarkerFaceColor','red','Color','red')
text(WPCA2model.loads{1,1}(qc2,1),WPCA2model.loads{1,1}(qc2,2),strtok(A(qc2,1),'.'))

hold on

s1cal = cellfun('isempty', strfind(Data.info(WPCA2.cal,2),'S'));
s2cal = find(s1cal == 0);
plot(WPCA2model.loads{1,1}(s2cal,1),WPCA2model.loads{1,1}(s2cal,2),'s','MarkerFaceColor','blue','Color','blue')
text(WPCA2model.loads{1,1}(s2cal,1),WPCA2model.loads{1,1}(s2cal,2),strtok(A(s2cal,1),'.'))

d1cal = cellfun('isempty', strfind(Data.info(WPCA2.cal,2),'D'));
d2cal = find(d1cal == 0);
plot(WPCA2model.loads{1,1}(d2cal,1),WPCA2model.loads{1,1}(d2cal,2),'s','MarkerFaceColor','green','Color','green')
text(WPCA2model.loads{1,1}(d2cal,1),WPCA2model.loads{1,1}(d2cal,2),strtok(A(d2cal,1),'.'))

B = Data.info(WPCA2.val,1);
s1val = cellfun('isempty', strfind(Data.info(WPCA2.val,2),'S'));
s2val = find(s1val == 0);
plot(WPCA2pred.loads{1,1}(s2val,1),WPCA2pred.loads{1,1}(s2val,2),'*','MarkerFaceColor','blue','Color','blue')
text(WPCA2pred.loads{1,1}(s2val,1),WPCA2pred.loads{1,1}(s2val,2),strtok(B(s2val,1),'.'))

d1val = cellfun('isempty', strfind(Data.info(WPCA2.val,2),'D'));
d2val = find(d1val == 0);
plot(WPCA2pred.loads{1,1}(d2val,1),WPCA2pred.loads{1,1}(d2val,2),'*','MarkerFaceColor','green','Color','green')
text(WPCA2pred.loads{1,1}(d2val,1),WPCA2pred.loads{1,1}(d2val,2),strtok(B(d2val,1),'.'))

legend('QC','S','D','Sval','Dval')



xlabel(['PC1 ','(',num2str(round(WPCA2model.detail.ssq(1,3),2)),'%',')'])
ylabel(['PC2 ','(',num2str(round(WPCA2model.detail.ssq(2,3),2)),'%',')'])

clear qc qc2 s1cal s2cal d1cal d2cal s1val s2val d1val d2val A B


%% EXTRACTING AND ORGANIZING THE LOADINGS from the PCA-W model (Pls-Toolbox)

%Creating the index
WPCA2B.ind = size(Xal{1,1},1)*size(Xal{1,1},2);
for i = 2:length(sic_groups);
    a = size(Xal{1,i},1)*size(Xal{1,i},2);
    WPCA2B.ind = [WPCA2B.ind WPCA2B.ind(i-1)+a];
end
clear a i

%% Extracting the LOADINGS

% Replacing NaN elements in WPCA.W

o = find(isnan(WPCA2.W) == 1);
WPCA2.W(1,o) = zeros(1,length(o));
clear o

W = reshape(WPCA2.W(1,1:size(Xal{1,1},1)*size(Xal{1,1},2)),size(Xal{1,1},1),size(Xal{1,1},2));

for j = 1:size(WPCA2model.loads{2,1},2);
    A = reshape(WPCA2model.loads{2,1}(1:size(Xal{1,1},1)*size(Xal{1,1},2),j),size(Xal{1,1},1),size(Xal{1,1},2));
    WPCA2.loads{1}(:,:,j) = A./W;
    for k = 1:size(WPCA2.loads{1}(:,:,j),1);
        o = find(isnan(WPCA2.loads{1}(k,:,j)) == 1);
        WPCA2.loads{1}(k,o,j) = zeros(1,length(o));
    end
    clear k o
end 
clear j A W

for i = 2:length(WPCA2.ind);
    W = reshape(WPCA2.W(1,WPCA2.ind(i-1)+1:WPCA2.ind(i)),size(Xal{1,i},1),size(Xal{1,i},2));
    
    for j = 1:size(WPCA2model.loads{2,1},2);
        A = reshape(WPCA2model.loads{2,1}(WPCA2.ind(i-1)+1:WPCA2.ind(i),j),size(Xal{1,i},1),size(Xal{1,i},2));
        WPCA2.loads{i}(:,:,j) = A./W;
        for k = 1:size(WPCA2.loads{i}(:,:,j),1);
            o = find(isnan(WPCA2.loads{i}(k,:,j)) == 1);
            WPCA2.loads{i}(k,o,j) = zeros(1,length(o));
        end
    clear k o
    end
    clear A W j
end
clear i

%% %%%%%%%%% CALCULATING THE CORRELATION LOADINGS %%%%%%%%%%%%%%%%%%%%%%%%

mc = Data.XalunDS.data(WPCA2.cal,:) - ones(size(Data.XalunDS.data(WPCA2.cal,:),1),1)*mean(Data.XalunDS.data(WPCA2.cal,:));
for i = 1:length(WPCA2.loads);
    mcsub = mc(:,WPCA2.ind(i)+1:WPCA2.ind(i+1));
    for j = 1:size(WPCA2.loads{i},3);
        load = reshape(WPCA2.loads{i}(:,:,j),1,size(WPCA2.loads{i}(:,:,j),1)*size(WPCA2.loads{i}(:,:,j),2));
        a = WPCA2model.loads{1,1}(:,j)'*WPCA2model.loads{1,1}(:,j);
        for k = 1:length(load);
            corrload(k) = load(1,k).*(sqrt(a)/sqrt(mcsub(:,k)'*mcsub(:,k)));
        end
        WPCA2.corrload{i}(:,:,j) = reshape(corrload,size(WPCA2.loads{i}(:,:,j),1),size(WPCA2.loads{i}(:,:,j),2));
    end 
    clear load a k corrload
end
clear mcsub mc i j


%% Ploting the loadings 
i = 8;
pc = 1;

figure;imagesc(Batch2.rt1{2}(rtLim{i}(1,1):rtLim{i}(1,2),1),Batch2.rt2{2}(rtLim{i}(2,1):rtLim{i}(2,2),1), WPCA2.loads{i}(:,:,pc))
axis xy

a = find(Batch2.acqT{2}./60 >= Batch2.rt1{2}(rtLim{i}(1,1)),1,'first');

b = Batch2.acqT{2}(a:a-1+size(WPCA2.loads{i},1)*size(WPCA2.loads{i},2))./60;

figure;plot(b,reshape(WPCA2.loads{i}(:,:,pc),1,size(WPCA2.loads{i},1)*size(WPCA2.loads{i},2))');
clear i pc a b


%% Defining Calibration and Validation Sets - WLS-PCA using MILES algorithm

WPCA.val = find(cell2mat(Data.info(:,3))); %Samples in Data.info(:,3) extracted for validation

WPCA.cal = find(cell2mat(Data.info(:,3)) == 0);
%%
% Computing WLS-PCA with MILES algorithm only in the calibration set

[WPCA2B.T,WPCA2B.P,WPCA2B.offset,WPCA2B.Loss,WPCA2B.Count]=milespca(Data.Xalun2(WPCA2B.cal,:),sparse(diag(WPCA2B.W)),3,1);


% Computing explained variances

for i = 1:size(WPCA1b.T,2);
        Dmc = Data.Xalun2(WPCA1b.cal,:)-ones(size(Data.Xalun2(WPCA1b.cal,:),1),1)*WPCA1b.offset;
        xcalc = WPCA1b.T(:,i)*WPCA1b.P(:,i)';
        var(i) = sum(sum(xcalc.^2))*100/sum(sum(Dmc.^2));
        WPCA1b.cumvar = cumsum(var);
end
clear xcalc var Dmc

figure;plot(WPCA2B.T(:,1),WPCA2B.T(:,2),'o')
text(WPCA2B.T(:,1),WPCA2B.T(:,2),strtok(Data.info(WPCA2B.cal,1),'.'))

% Calculating scores for the validation samples
WPCA2B.Tval = (Data.Xalun2(WPCA2B.val,:) - ones(size(Data.Xalun2(WPCA2B.val,:),1),1)*WPCA2B.offset)*diag(WPCA2B.W)*WPCA2B.P*inv(WPCA2B.P'*diag(WPCA2B.W)*WPCA2B.P);


%% Plotting the PC1 x PC2 scores (no validation)

qc = cellfun('isempty', strfind(Data.info(:,2),'QC'));
qc2 = find(qc == 0);
A = Data.info(:,1);
figure;plot(WPCA2B.T(qc2,1),WPCA2B.T(qc2,2),'o','MarkerFaceColor','red','Color','red')
text(WPCA2B.T(qc2,1),WPCA2B.T(qc2,2),strtok(A(qc2,1),'.'))

hold on

s1 = cellfun('isempty', strfind(Data.info(:,2),'S'));
s2 = find(s1 == 0);
plot(WPCA2B.T(s2,1),WPCA2B.T(s2,2),'s','MarkerFaceColor','blue','Color','blue')
text(WPCA2B.T(s2,1),WPCA2B.T(s2,2),strtok(A(s2,1),'.'))

d1 = cellfun('isempty', strfind(Data.info(:,2),'D'));
d2 = find(d1 == 0);
plot(WPCA2B.T(d2,1),WPCA2B.T(d2,2),'s','MarkerFaceColor','green','Color','green')
text(WPCA2B.T(d2,1),WPCA2B.T(d2,2),strtok(A(d2,1),'.'))

legend('QC','S','D')

xlabel('PC1')
ylabel('PC2')

clear qc qc2 s1 s2 d1 d2 A B

%% Plotting the PC1 x PC2 scores (with validation) - LOCAL PCA
i = 3;

qc = cellfun('isempty', strfind(Data.info(LocalPCA.cal,2),'QC'));
qc2 = find(qc == 0);
A = Data.info(LocalPCA.cal,1);
figure;plot(LocalPCA.T{i}(qc2,1),LocalPCA.T{i}(qc2,2),'o','MarkerFaceColor','red','Color','red')
text(LocalPCA.T{i}(qc2,1),LocalPCA.T{i}(qc2,2),strtok(A(qc2,1),'.'))

hold on

s1cal = cellfun('isempty', strfind(Data.info(LocalPCA.cal,2),'S'));
s2cal = find(s1cal == 0);
plot(LocalPCA.T{i}(s2cal,1),LocalPCA.T{i}(s2cal,2),'s','MarkerFaceColor','blue','Color','blue')
text(LocalPCA.T{i}(s2cal,1),LocalPCA.T{i}(s2cal,2),strtok(A(s2cal,1),'.'))

d1cal = cellfun('isempty', strfind(Data.info(LocalPCA.cal,2),'D'));
d2cal = find(d1cal == 0);
plot(LocalPCA.T{i}(d2cal,1),LocalPCA.T{i}(d2cal,2),'s','MarkerFaceColor','green','Color','green')
text(LocalPCA.T{i}(d2cal,1),LocalPCA.T{i}(d2cal,2),strtok(A(d2cal,1),'.'))

B = Data.info(LocalPCA.val,1);
s1val = cellfun('isempty', strfind(Data.info(LocalPCA.val,2),'S'));
s2val = find(s1val == 0);
plot(LocalPCA.Tpred{i}(s2val,1),LocalPCA.Tpred{i}(s2val,2),'*','MarkerFaceColor','blue','Color','blue')
text(LocalPCA.Tpred{i}(s2val,1),LocalPCA.Tpred{i}(s2val,2),strtok(B(s2val,1),'.'))

d1val = cellfun('isempty', strfind(Data.info(LocalPCA.val,2),'D'));
d2val = find(d1val == 0);
plot(LocalPCA.Tpred{i}(d2val,1),LocalPCA.Tpred{i}(d2val,2),'*','MarkerFaceColor','green','Color','green')
text(LocalPCA.Tpred{i}(d2val,1),LocalPCA.Tpred{i}(d2val,2),strtok(B(d2val,1),'.'))

legend('QC','S','D','Sval','Dval')

xlabel(['PC1 ','(',num2str(round(LocalPCA.SSQ{i}(1,3),2)),'%',')'])
ylabel(['PC2 ','(',num2str(round(LocalPCA.SSQ{i}(2,3),2)),'%',')'])
title(['SIC = ',num2str(sic_groups{i})])
clear qc qc2 s1cal s2cal d1cal d2cal s1val s2val d1val d2val A B



%% Plotting the PC1 x PC2 scores (no validation)

qc = cellfun('isempty', strfind(Data.info(:,2),'QC'));
qc2 = find(qc == 0);
A = Data.info(:,1);
figure;plot(WPCA2B.T(qc2,1),WPCA2B.T(qc2,2),'o','MarkerFaceColor','red','Color','red')
text(WPCA2B.T(qc2,1),WPCA2B.T(qc2,2),strtok(A(qc2,1),'.'))

hold on

s1 = cellfun('isempty', strfind(Data.info(:,2),'S'));
s2 = find(s1 == 0);
plot(WPCA2B.T(s2,1),WPCA2B.T(s2,2),'s','MarkerFaceColor','blue','Color','blue')
text(WPCA2B.T(s2,1),WPCA2B.T(s2,2),strtok(A(s2,1),'.'))

d1 = cellfun('isempty', strfind(Data.info(:,2),'D'));
d2 = find(d1 == 0);
plot(WPCA2B.T(d2,1),WPCA2B.T(d2,2),'s','MarkerFaceColor','green','Color','green')
text(WPCA2B.T(d2,1),WPCA2B.T(d2,2),strtok(A(d2,1),'.'))

legend('QC','S','D')

xlabel('PC1')
ylabel('PC2')

clear qc qc2 s1 s2 d1 d2 A B


%% INTERPRETING THE LOADINGS FROM PCA-W (W = RSD^-1)

%% Sorting the loadings according to relative importance

%% Establishing the threshold to sort the loadings
%threshold1 = intensity to distiguish the peaks from noise
%threshold2 = overall mean of the loadings in the i-th PC

pc = 2;
overall_loads.threshold1(pc) = 5e-5;
loadcor = abs(WPCA2model.loads{2,1}(:,pc));
loadcor = loadcor./WPCA2.W';

c = find(isnan(loadcor) == 1);
loadcor(c,1) = zeros(length(c),1);

a = find(loadcor > overall_loads.threshold1(pc));
a2 = abs(loadcor(a,1));
overall_loads.threshold2(pc) = mean(a2);
findpeaks(abs(loadcor)','MinPeakHeight',overall_loads.threshold2(pc));
clear a a2

clear c i j pc loadcor

%% Computing contribution of the loadings
pc = 2;
for i = 1:length(sic_groups);
    b = abs(reshape(WPCA2.loads{i}(:,:,pc),1,size(WPCA2.loads{i},1)*size(WPCA2.loads{i},2)));
    [int,~,~,~] = findpeaks(b,'MinPeakHeight',overall_loads.threshold2(pc));
    [int2,~,~,~] = findpeaks(b,'MinPeakHeight',overall_loads.threshold1(pc));
    overall_loads.peaksratio(pc,i) = length(int)/length(int2);
    clear b int int2
end


%%
figure; for i = 1:size(overall_loads.peaksratio,1);   
bar(overall_loads.peaksratio(i,:)),hold on,axis tight,end
%%
for j = 1:length(overall_loads.peaksratio);
    if isnan(overall_loads.peaksratio(j)) == 0;
        figure;plot(abs(reshape(WPCA2.loads{j}(:,:,pc),1,size(WPCA2.loads{j},1)*size(WPCA2.loads{j},2))));
        title(['SIC=',num2str(sic_groups{j})]);
    end
end


%% %%%%%%%%%%%%%%% BUILDING LOCAL AND SUCESSIVE W-PCA MODELS %%%%%%%%%%%%%
%%%%%%%%%%   EVALUATION OF VARIABLES MORE SUSCEPTIBLE TO WEATHERING %%%%%%
%%%%%%%%%  For the SIC from each PAHs, only a global model is computed %%%


%% Running the sucessive W-PCA models

% Fitting the overall W-PCA models in each SIC

h = waitbar(0,'Computing individual SICs W-PCA models...');

for i = 1:length(Data.ind)-1;
    waitbar(i / length(Data.ind))
    
    X = Data.Xalun2N1(:,Data.ind(i)+1:Data.ind(i+1));
    w = LocalPCA.W(1,Data.ind(i)+1:Data.ind(i+1));

    % Fitting the model

    Xcal = X(LocalPCA.cal,:);
    Xval = X(LocalPCA.val,:);

    optcv = crossval('options');
    optcv.display = 'off';
    optcv.plots = 'none';
    optcv.rmsec = 'no';
    optcv.structureoutput = 'no';
    optcv.pcacvi = {'vet' 7};
    optcv.preprocessing = 1; %mean center

    [~,~,rmsecv,~,~,~,~] = crossval(Xcal,[],'pca',{'vet' 7},5,optcv);
    a = find(rmsecv == min(rmsecv),1,'first');
    LocalPCA.rmsecv{i} = rmsecv(a); clear rmsecv

    optpca = pca('options');
    optpca.display = 'off';
    optpca.plots = 'none';
    optpca.preprocessing = {preprocess('default','mean center') []};

    model = pca(Xcal,a,optpca);
    LocalPCA.T{i} = model.loads{1,1};
    for k = 1:size(model.loads{2,1},2);
        LocalPCA.P{i}(:,k) = model.loads{2,1}(:,k)./w';
        c = find(isnan(LocalPCA.P{i}(:,k)) == 1);
        LocalPCA.P{i}(c,k) = zeros(length(c),1);
    end
  
    LocalPCA.SSQ{i} = model.detail.ssq(1:a,:);

    pred = pca(Xval,model,optpca);
    LocalPCA.Tpred{i} = pred.loads{1,1};
    clear model pred a X Xcal Xval w k c
end
close(h)
clear h i

%% %%%%%% COMPUTING STATISTICS FROM THE SCORES WITHIN EACH SIC %%%%%%%%%%%
%         Ranking the importance of each SIC for the PCA model


% Computing the averaged variances of the mean PCA Scores and the Scores from
% replicates


for i = 1:length(LocalPCA.T);
    
    LocalPCA.Tmean{i}(1,:) = mean(LocalPCA.T{i}([1,2,15,26,37],:)); %Averaged QCs
    LocalPCA.Tmean{i}(2,:) = mean([LocalPCA.T{i}(3,:);LocalPCA.Tpred{i}(6,:)]);
    LocalPCA.Tmean{i}(3:7,:) = LocalPCA.T{i}(4:8,:);
    LocalPCA.Tmean{i}(8,:) = mean([LocalPCA.T{i}(9,:);LocalPCA.Tpred{i}(1,:)]);
    LocalPCA.Tmean{i}(9:11,:) = LocalPCA.T{i}(10:12,:);
    LocalPCA.Tmean{i}(12,:) = mean([LocalPCA.T{i}(13,:);LocalPCA.Tpred{i}(2,:)]);
    LocalPCA.Tmean{i}(13:15,:) = LocalPCA.T{i}([14,16,17],:);
    LocalPCA.Tmean{i}(16,:) = mean([LocalPCA.T{i}(18,:);LocalPCA.Tpred{i}(3,:)]);
    LocalPCA.Tmean{i}(17,:) = LocalPCA.T{i}(19,:);
    LocalPCA.Tmean{i}(18,:) = mean([LocalPCA.T{i}(20,:);LocalPCA.Tpred{i}(4,:)]);
    LocalPCA.Tmean{i}(19:20,:) = LocalPCA.T{i}(21:22,:);
    LocalPCA.Tmean{i}(21,:) = mean([LocalPCA.T{i}(23,:);LocalPCA.Tpred{i}(5,:)]);
    LocalPCA.Tmean{i}(22:24,:) = LocalPCA.T{i}([24,25,27],:);
    LocalPCA.Tmean{i}(25,:) = mean([LocalPCA.T{i}(28,:);LocalPCA.Tpred{i}(7,:)]);
    LocalPCA.Tmean{i}(26,:) = mean([LocalPCA.T{i}(29,:);LocalPCA.Tpred{i}(8,:)]);
    LocalPCA.Tmean{i}(27:28,:) = LocalPCA.T{i}(30:31,:);
    LocalPCA.Tmean{i}(29,:) = mean([LocalPCA.T{i}(32,:);LocalPCA.Tpred{i}(9,:)]);
    LocalPCA.Tmean{i}(30:31,:) = LocalPCA.T{i}(33:34,:);
    LocalPCA.Tmean{i}(32,:) = mean([LocalPCA.T{i}(35,:);LocalPCA.Tpred{i}(10,:)]);
    LocalPCA.Tmean{i}(33,:) = LocalPCA.T{i}(36,:);
    
    LocalPCA.DiagSIC{i}(:,1) = mean(var(LocalPCA.Tmean{i})'.*(LocalPCA.SSQ{i}(:,3)./100));
    
end

for i = 1:length(LocalPCA.T);
    X1 = mean(var(LocalPCA.T{i}([1,2,15,26,37],:))'.*(LocalPCA.SSQ{i}(:,3)./100));
    X2 = mean(var([LocalPCA.T{i}(3,:);LocalPCA.Tpred{i}(6,:)])'.*(LocalPCA.SSQ{i}(:,3)./100));
    X3 = mean(var([LocalPCA.T{i}(9,:);LocalPCA.Tpred{i}(1,:)])'.*(LocalPCA.SSQ{i}(:,3)./100));
    X4 = mean(var([LocalPCA.T{i}(13,:);LocalPCA.Tpred{i}(2,:)])'.*(LocalPCA.SSQ{i}(:,3)./100));
    X5 = mean(var([LocalPCA.T{i}(18,:);LocalPCA.Tpred{i}(3,:)])'.*(LocalPCA.SSQ{i}(:,3)./100));
    X6 = mean(var([LocalPCA.T{i}(20,:);LocalPCA.Tpred{i}(4,:)])'.*(LocalPCA.SSQ{i}(:,3)./100));
    X7 = mean(var([LocalPCA.T{i}(23,:);LocalPCA.Tpred{i}(5,:)])'.*(LocalPCA.SSQ{i}(:,3)./100));
    X8 = mean(var([LocalPCA.T{i}(28,:);LocalPCA.Tpred{i}(7,:)])'.*(LocalPCA.SSQ{i}(:,3)./100));;
    X9 = mean(var([LocalPCA.T{i}(29,:);LocalPCA.Tpred{i}(8,:)])'.*(LocalPCA.SSQ{i}(:,3)./100));
    X10 = mean(var([LocalPCA.T{i}(32,:);LocalPCA.Tpred{i}(9,:)])'.*(LocalPCA.SSQ{i}(:,3)./100));
    X11 = mean(var([LocalPCA.T{i}(35,:);LocalPCA.Tpred{i}(10,:)])'.*(LocalPCA.SSQ{i}(:,3)./100));;

    
    LocalPCA.DiagSIC{i}(:,2) = ((4*X1)+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11)/14;
    
    clear X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11

end

%% Computing F-ratio Statistitics between Model vs Replicates variances

for i = 1:length(LocalPCA.T);
    
    LocalPCA.DiagSIC{i}(:,3) = LocalPCA.DiagSIC{i}(:,1)./LocalPCA.DiagSIC{i}(:,2);
    
    LocalPCA.DiagSIC{i}(:,4) = finv(0.95,size(LocalPCA.Tmean{i},1),14);
    
end
%%
for i = 1:length(LocalPCA.T);
   a(:,i) = [LocalPCA.DiagSIC{i}(:,3);LocalPCA.DiagSIC{i}(:,4)];
end
figure;plot(a(1,:)'),hold on, plot(a(2,:)'),clear a i

%% PERFORMING W-PCA ONLY IN THE SIGNIFICANT SIGNALS

%% Selecting the SICs
WPCA.sics = [1:15,17:23,25:32];

X = [];

for j = 1:length(WPCA.sics);

    X = [X Data.Xalun2N1(:,Data.ind(WPCA.sics(j))+1:Data.ind(WPCA.sics(j)+1))];
    
end

clear j
%% Running the WPCA model

WPCA.cal = LocalPCA.cal;
WPCA.val = LocalPCA.val;

% Creating Weights

X1 = X([2,1,17,32,47],:);   
w1 = 1./(std(X1)./mean(X1));
c = find(isnan(w1) == 1);
w1(1,c) = zeros(1,length(c));

WPCA.W = w1;
clear X1 c w1


for k = 1:size(X,1);
    X2(k,:) = X(k,:).*WPCA.W; %downscaling variables
end

Xcal = X2(WPCA.cal,:);
Xval = X2(WPCA.val,:);

[~,~,rmsecv,~,~,~,~] = crossval(Xcal,[],'pca',{'vet' 7},4,optcv);
a = find(rmsecv == min(rmsecv),1,'first');
clear rmsecv

WPCAmodel = pca(Xcal,a,optpca);
for k = 1:size(WPCAmodel.loads{2,1},2);
    WPCA.loads(:,k) = WPCAmodel.loads{2,1}(:,k)./WPCA.W'; %upscaling back for interpretation
    c = find(isnan(WPCA.loads(:,k)) == 1);
    WPCA.loads(c,k) = zeros(length(c),1);
    c = find(isfinite(WPCA.loads(:,k)) == 0);
    WPCA.loads(c,k) = zeros(length(c),1);
end
  
    WPCApred = pca(Xval,WPCAmodel,optpca);
    clear a Xcal Xval k c i X X2

    
 %% Plotting the WPC1 x WPC2 scores (with validation)

qc = cellfun('isempty', strfind(Data.info(WPCA.cal,2),'QC'));
qc2 = find(qc == 0);
A = Data.info(WPCA.cal,1);
figure;plot(WPCAmodel.loads{1,1}(qc2,1),WPCAmodel.loads{1,1}(qc2,2),'o','MarkerFaceColor','red','Color','red')
text(WPCAmodel.loads{1,1}(qc2,1),WPCAmodel.loads{1,1}(qc2,2),strtok(A(qc2,1),'.'))

hold on

s1cal = cellfun('isempty', strfind(Data.info(WPCA.cal,2),'S'));
s2cal = find(s1cal == 0);
plot(WPCAmodel.loads{1,1}(s2cal,1),WPCAmodel.loads{1,1}(s2cal,2),'s','MarkerFaceColor','blue','Color','blue')
text(WPCAmodel.loads{1,1}(s2cal,1),WPCAmodel.loads{1,1}(s2cal,2),strtok(A(s2cal,1),'.'))

d1cal = cellfun('isempty', strfind(Data.info(WPCA.cal,2),'D'));
d2cal = find(d1cal == 0);
plot(WPCAmodel.loads{1,1}(d2cal,1),WPCAmodel.loads{1,1}(d2cal,2),'s','MarkerFaceColor','green','Color','green')
text(WPCAmodel.loads{1,1}(d2cal,1),WPCAmodel.loads{1,1}(d2cal,2),strtok(A(d2cal,1),'.'))

B = Data.info(WPCA.val,1);
s1val = cellfun('isempty', strfind(Data.info(WPCA.val,2),'S'));
s2val = find(s1val == 0);
plot(WPCApred.loads{1,1}(s2val,1),WPCApred.loads{1,1}(s2val,2),'*','MarkerFaceColor','blue','Color','blue')
text(WPCApred.loads{1,1}(s2val,1),WPCApred.loads{1,1}(s2val,2),strtok(B(s2val,1),'.'))

d1val = cellfun('isempty', strfind(Data.info(WPCA.val,2),'D'));
d2val = find(d1val == 0);
plot(WPCApred.loads{1,1}(d2val,1),WPCApred.loads{1,1}(d2val,2),'*','MarkerFaceColor','green','Color','green')
text(WPCApred.loads{1,1}(d2val,1),WPCApred.loads{1,1}(d2val,2),strtok(B(d2val,1),'.'))

legend('QC','S','D','Sval','Dval')



xlabel(['PC1 ','(',num2str(round(WPCAmodel.detail.ssq(1,3),2)),'%',')'])
ylabel(['PC2 ','(',num2str(round(WPCAmodel.detail.ssq(2,3),2)),'%',')'])

clear qc qc2 s1cal s2cal d1cal d2cal s1val s2val d1val d2val A B



%% Computing confidence intervals for the Scores from the QCs errors

    WPCA.Tmean(1,:) = mean(WPCAmodel.loads{1,1}([1,2,15,26,37],:)); %Averaged QCs
    WPCA.Tmean(2,:) = mean([WPCAmodel.loads{1,1}(3,:);WPCApred.loads{1,1}(6,:)]);
    WPCA.Tmean(3:7,:) = WPCAmodel.loads{1,1}(4:8,:);
    WPCA.Tmean(8,:) = mean([WPCAmodel.loads{1,1}(9,:);WPCApred.loads{1,1}(1,:)]);
    WPCA.Tmean(9:11,:) = WPCAmodel.loads{1,1}(10:12,:);
    WPCA.Tmean(12,:) = mean([WPCAmodel.loads{1,1}(13,:);WPCApred.loads{1,1}(2,:)]);
    WPCA.Tmean(13:15,:) = WPCAmodel.loads{1,1}([14,16,17],:);
    WPCA.Tmean(16,:) = mean([WPCAmodel.loads{1,1}(18,:);WPCApred.loads{1,1}(3,:)]);
    WPCA.Tmean(17,:) = WPCAmodel.loads{1,1}(19,:);
    WPCA.Tmean(18,:) = mean([WPCAmodel.loads{1,1}(20,:);WPCApred.loads{1,1}(4,:)]);
    WPCA.Tmean(19:20,:) = WPCAmodel.loads{1,1}(21:22,:);
    WPCA.Tmean(21,:) = mean([WPCAmodel.loads{1,1}(23,:);WPCApred.loads{1,1}(5,:)]);
    WPCA.Tmean(22:24,:) = WPCAmodel.loads{1,1}([24,25,27],:);
    WPCA.Tmean(25,:) = mean([WPCAmodel.loads{1,1}(28,:);WPCApred.loads{1,1}(7,:)]);
    WPCA.Tmean(26,:) = mean([WPCAmodel.loads{1,1}(29,:);WPCApred.loads{1,1}(8,:)]);
    WPCA.Tmean(27:28,:) = WPCAmodel.loads{1,1}(30:31,:);
    WPCA.Tmean(29,:) = mean([WPCAmodel.loads{1,1}(32,:);WPCApred.loads{1,1}(9,:)]);
    WPCA.Tmean(30:31,:) = WPCAmodel.loads{1,1}(33:34,:);
    WPCA.Tmean(32,:) = mean([WPCAmodel.loads{1,1}(35,:);WPCApred.loads{1,1}(10,:)]);
    WPCA.Tmean(33,:) = WPCAmodel.loads{1,1}(36,:);
    
    % Computing variances of each PC
    
    X1 = var(WPCAmodel.loads{1,1}([1,2,15,26,37],:));
    X2 = var([WPCAmodel.loads{1,1}(3,:);WPCApred.loads{1,1}(6,:)]);
    X3 = var([WPCAmodel.loads{1,1}(9,:);WPCApred.loads{1,1}(1,:)]);
    X4 = var([WPCAmodel.loads{1,1}(13,:);WPCApred.loads{1,1}(2,:)]);
    X5 = var([WPCAmodel.loads{1,1}(18,:);WPCApred.loads{1,1}(3,:)]);
    X6 = var([WPCAmodel.loads{1,1}(20,:);WPCApred.loads{1,1}(4,:)]);
    X7 = var([WPCAmodel.loads{1,1}(23,:);WPCApred.loads{1,1}(5,:)]);
    X8 = var([WPCAmodel.loads{1,1}(28,:);WPCApred.loads{1,1}(7,:)]);
    X9 = var([WPCAmodel.loads{1,1}(29,:);WPCApred.loads{1,1}(8,:)]);
    X10 = var([WPCAmodel.loads{1,1}(32,:);WPCApred.loads{1,1}(9,:)]);
    X11 = var([WPCAmodel.loads{1,1}(35,:);WPCApred.loads{1,1}(10,:)]);

    for i = 1:size(X1,2);
        
        WPCA.stdpooled(i) = sqrt(((4*X1(1,i))+X2(1,i)+X3(1,i)+X4(1,i)+X5(1,i)+...
            X6(1,i)+X7(1,i)+X8(1,i)+X9(1,i)+X10(1,i)+X11(1,i))/14);
        
    end
    
        t = tinv(0.95,14); % t-value for 14 dfs at 95% conf. level

    
    for i = 1:size(X1,2);
    
        for j = 1:size(WPCA.Tmean,1);
            WPCA.Tconfint(j,i) = abs((WPCA.Tmean(j,i)+(WPCA.stdpooled(i)*t/sqrt(25))) -...
                (WPCA.Tmean(j,i)-(WPCA.stdpooled(i)*t/sqrt(25)))); % Computing confidence intervals
        end
    end 
      
    clear X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 i j
    
%% Plotting mean Scores and Confidence Intervals
    
samples = Data.info([1,3:9,11:14,16,18:20,22:23,25:27,29:30,33:34,36,38:40,42:44,46],1:2);

s1cal = cellfun('isempty', strfind(samples(:,2),'S'));
s2cal = find(s1cal == 0);
figure;errorbar(WPCA.Tmean(s2cal,1),WPCA.Tmean(s2cal,2),WPCA.Tconfint(s2cal,2),'s','MarkerFaceColor','blue','Color','blue')
hold on
herrorbar(WPCA.Tmean(s2cal,1),WPCA.Tmean(s2cal,2),WPCA.Tconfint(s2cal,1),'blue')
text(WPCA.Tmean(s2cal,1),WPCA.Tmean(s2cal,2),strtok(samples(s2cal,1),'.'))

d1cal = cellfun('isempty', strfind(samples(:,2),'D'));
d2cal = find(d1cal == 0);
errorbar(WPCA.Tmean(d2cal,1),WPCA.Tmean(d2cal,2),WPCA.Tconfint(d2cal,2),'s','MarkerFaceColor','green','Color','green')
herrorbar(WPCA.Tmean(d2cal,1),WPCA.Tmean(d2cal,2),WPCA.Tconfint(d2cal,1),'green')
text(WPCA.Tmean(d2cal,1),WPCA.Tmean(d2cal,2),strtok(samples(d2cal,1),'.'))


%legend('QC','S','D')

xlabel(['PC1 ','(',num2str(round(WPCAmodel.detail.ssq(1,3),2)),'%',')'])
ylabel(['PC2 ','(',num2str(round(WPCAmodel.detail.ssq(2,3),2)),'%',')'])

clear s1cal s2cal d1cal d2cal samples


%% Calculating distance amongst Scores

WPCA2.Cases2.c1.samples = {'3182-2';'3211-5';'3243-12';'3243-7'};
WPCA2.Cases2.c2.samples = {'3188-2';'3211-4';'3243-8'};
WPCA2.Cases2.c3.samples = {'3184-2';'3211-4';'3243-8'};
WPCA2.Cases2.c4.samples = {'3184-1';'3211-4';'3243-8'};

WPCA2.Cases2.c1.Tmean = WPCA2.Tmean([12,23,32,15],:);
WPCA2.Cases2.c2.Tmean = WPCA2.Tmean([6,2,14],:);
WPCA2.Cases2.c3.Tmean = WPCA2.Tmean([7,2,14],:);
WPCA2.Cases2.c4.Tmean = WPCA2.Tmean([8,2,14],:);


%% Computing probability of potential sources
for i = 1:4; %Considering only the first 5 PCs (more relevant)
    k = 1; %spill sample in WPCA2.Cases1.cn.samples
    
    for j = 2; %sources
        
        q = abs(WPCA2.Cases1.c1.Tmean(k,i) - WPCA2.Cases1.c1.Tmean(j,i))/(WPCA2.stdpooled(i)*sqrt(2));
        
        WPCA2.Cases1.c1.prob{j}(k,i) = 2*(1 - tcdf(q,14));
        clear q
    end
    
end
for j = 2;
    WPCA2.Cases1.c4.combprob(k,j) = sum(WPCA2.Cases1.c4.prob{j}(k,1:5)'.*(WPCA2model.detail.ssq(1:5,3)./100))/(WPCA2model.detail.ssq(5,4)/100);
end
clear i j k




%% PLOTTING LOADINGS

ind2 = size(Xal{1,1},1)*size(Xal{1,1},2);
for i = 2:length(WPCA.sics);
    a = size(Xal{1,WPCA.sics(i)},1)*size(Xal{1,WPCA.sics(i)},2);
    ind2 = [ind2 ind2(i-1)+a];
end
ind2 = [0 ind2];
clear a i

sicpos = [16];
X = zeros(size(Batch2.SIC_GPS{2},1),size(Batch2.SIC_GPS{2},1));

pc = 2;

for i = 1:length(sicpos);
    X2 = WPCA.loads(ind2(sicpos(i))+1:ind2(sicpos(i)+1),pc);
    X2 = reshape(X2,size(Xal{1,WPCA.sics(sicpos(i))},1),size(Xal{1,WPCA.sics(sicpos(i))},2));

    X(rtLim{WPCA.sics(sicpos(i))}(2,1):rtLim{WPCA.sics(sicpos(i))}(2,2),rtLim{WPCA.sics(sicpos(i))}(1,1):rtLim{WPCA.sics(sicpos(i))}(1,2)) = X2;
end

figure;imagesc(X)
axis xy

clear ind2 sicpos X X2 i pc
  
%% PERFORMING OIL-SOURCE CORRELATION WHILE HANDLING EVAPORATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computing SUCESSIVE WPCA models CONSIDERING ALL SAMPLES.
% W computed only from the QCs
% Euclidian norm of the entire row-vector


%% Selecting the signifcant SICs (from the LocalPCAs)

%% Defining the variable intervals for local modelling

% Intervals based on alkanes rt in the aligned 2D GC chromatograms

WPCA1_eva.sics = WPCA.sics;
WPCA1_eva.rt1int = [9.9,29.94,32.63,34.85,37.18,39.17,41.03,43.02,...
    44.88,46.52,48.15];

%% Running SUCESSIVE WPCA models eliminating evaporation

% CALCULATING SUCESSIVE WPCA MODELS

%WPCA1_eva.val = WPCA.val;
%WPCA1_eva.cal = WPCA.cal;


for i = 2:length(WPCA1_eva.rt1int);
    
    Xmodel = [];
    w = [];
    
    for j = 1:length(WPCA1_eva.sics);

    a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(j)}(1,1):rtLim{WPCA1_eva.sics(j)}(1,2),1) >= WPCA1_eva.rt1int(i),1,'first');

        for k = 1:size(Xal3,1);
            A = Xal3{k,j}(:,a:end);
            X(k,:) = reshape(A,1,size(A,1)*size(A,2));
        end
   
    Xmodel = [Xmodel X];
    
    clear X A
    end  
    
    Xmodel = normr(Xmodel); % Normalizing entire row to euclidian norm = 1
    
    
    %Extracting the corresponding weights
   
    X1 = Xmodel([2,1,17,32,47],:);
    
    w = 1./(std(X1)./mean(X1));

    c = find(isnan(w) == 1);
    w(1,c) = zeros(1,length(c));

    clear X1 c
        
        
    % Calculating the models
    
    for k = 1:size(Xmodel,1);
        X(k,:) = Xmodel(k,:).*w; %downscaling variables
    end
    
    Xcal = X(WPCA1_eva.cal,:);
    Xval = X(WPCA1_eva.val,:);
    
    [~,~,rmsecv,~,~,~,~] = crossval(Xcal,[],'pca',{'vet' 5},3,optcv);
    a = find(rmsecv == min(rmsecv),1,'first');
    
    model = pca(Xcal,a,optpca);
    
    WPCA1_eva.models{i}.T = model.loads{1,1}; % extracting scores
    
    for k = 1:size(model.loads{2,1},2);
        WPCA1_eva.models{i}.P(:,k) = model.loads{2,1}(:,k)./w';
        c = find(isnan(WPCA1_eva.models{i}.P(:,k)) == 1);
        WPCA1_eva.models{i}.P(c,k) = zeros(length(c),1);
        c = find(isfinite(WPCA1_eva.models{i}.P(:,k)) == 0);
        WPCA1_eva.models{i}.P(c,k) = zeros(length(c),1);
    end
  
    WPCA1_eva.models{i}.SSQ = model.detail.ssq(1:a,:);
    
    pred = pca(Xval,model,optpca);
    WPCA1_eva.models{i}.Tpred = pred.loads{1,1}; % extracting validation scores
    
    clear model pred k Xcal Xval a rmsecv Xmodel w X
    
end

clear a i j

 %% Plotting the PC1 x PC2 scores (with validation)
i = 5;

qc = cellfun('isempty', strfind(Data.info(WPCA1_eva.cal,2),'QC'));
qc2 = find(qc == 0);
A = Data.info(WPCA1_eva.cal,1);
figure;plot(WPCA1_eva.models{i}.T(qc2,1),WPCA1_eva.models{i}.T(qc2,2),'o','MarkerFaceColor','red','Color','red')
text(WPCA1_eva.models{i}.T(qc2,1),WPCA1_eva.models{i}.T(qc2,2),strtok(A(qc2,1),'.'))

hold on

s1cal = cellfun('isempty', strfind(Data.info(WPCA1_eva.cal,2),'S'));
s2cal = find(s1cal == 0);
plot(WPCA1_eva.models{i}.T(s2cal,1),WPCA1_eva.models{i}.T(s2cal,2),'s','MarkerFaceColor','blue','Color','blue')
text(WPCA1_eva.models{i}.T(s2cal,1),WPCA1_eva.models{i}.T(s2cal,2),strtok(A(s2cal,1),'.'))

d1cal = cellfun('isempty', strfind(Data.info(WPCA1_eva.cal,2),'D'));
d2cal = find(d1cal == 0);
plot(WPCA1_eva.models{i}.T(d2cal,1),WPCA1_eva.models{i}.T(d2cal,2),'s','MarkerFaceColor','green','Color','green')
text(WPCA1_eva.models{i}.T(d2cal,1),WPCA1_eva.models{i}.T(d2cal,2),strtok(A(d2cal,1),'.'))

B = Data.info(WPCA1_eva.val,1);
s1val = cellfun('isempty', strfind(Data.info(WPCA1_eva.val,2),'S'));
s2val = find(s1val == 0);
plot(WPCA1_eva.models{i}.Tpred(s2val,1),WPCA1_eva.models{i}.Tpred(s2val,2),'*','MarkerFaceColor','blue','Color','blue')
text(WPCA1_eva.models{i}.Tpred(s2val,1),WPCA1_eva.models{i}.Tpred(s2val,2),strtok(B(s2val,1),'.'))

d1val = cellfun('isempty', strfind(Data.info(WPCA1_eva.val,2),'D'));
d2val = find(d1val == 0);
plot(WPCA1_eva.models{i}.Tpred(d2val,1),WPCA1_eva.models{i}.Tpred(d2val,2),'*','MarkerFaceColor','green','Color','green')
text(WPCA1_eva.models{i}.Tpred(d2val,1),WPCA1_eva.models{i}.Tpred(d2val,2),strtok(B(d2val,1),'.'))

legend('QC','S','D','Sval','Dval')

xlabel(['PC1 ','(',num2str(round(WPCA1_eva.models{i}.SSQ(1,3),2)),'%',')'])
ylabel(['PC2 ','(',num2str(round(WPCA1_eva.models{i}.SSQ(2,3),2)),'%',')'])
title(['rt1D (min) > ',num2str(WPCA1_eva.rt1int(i))])
clear qc qc2 s1cal s2cal d1cal d2cal s1val s2val d1val d2val A B

    

%% Computing confidence intervals for the Scores from the QCs errors in the
% WPCA2_eva models

for i = 2:length(WPCA1_eva.models);

    WPCA1_eva.models{i}.Tmean(1,:) = mean(WPCA1_eva.models{i}.T([1,2,15,26,37],:)); %Averaged QCs
    WPCA1_eva.models{i}.Tmean(2,:) = mean([WPCA1_eva.models{i}.T(3,:);WPCA1_eva.models{i}.Tpred(6,:)]);
    WPCA1_eva.models{i}.Tmean(3:7,:) = WPCA1_eva.models{i}.T(4:8,:);
    WPCA1_eva.models{i}.Tmean(8,:) = mean([WPCA1_eva.models{i}.T(9,:);WPCA1_eva.models{i}.Tpred(1,:)]);
    WPCA1_eva.models{i}.Tmean(9:11,:) = WPCA1_eva.models{i}.T(10:12,:);
    WPCA1_eva.models{i}.Tmean(12,:) = mean([WPCA1_eva.models{i}.T(13,:);WPCA1_eva.models{i}.Tpred(2,:)]);
    WPCA1_eva.models{i}.Tmean(13:15,:) = WPCA1_eva.models{i}.T([14,16,17],:);
    WPCA1_eva.models{i}.Tmean(16,:) = mean([WPCA1_eva.models{i}.T(18,:);WPCA1_eva.models{i}.Tpred(3,:)]);
    WPCA1_eva.models{i}.Tmean(17,:) = WPCA1_eva.models{i}.T(19,:);
    WPCA1_eva.models{i}.Tmean(18,:) = mean([WPCA1_eva.models{i}.T(20,:);WPCA1_eva.models{i}.Tpred(4,:)]);
    WPCA1_eva.models{i}.Tmean(19:20,:) = WPCA1_eva.models{i}.T(21:22,:);
    WPCA1_eva.models{i}.Tmean(21,:) = mean([WPCA1_eva.models{i}.T(23,:);WPCA1_eva.models{i}.Tpred(5,:)]);
    WPCA1_eva.models{i}.Tmean(22:24,:) = WPCA1_eva.models{i}.T([24,25,27],:);
    WPCA1_eva.models{i}.Tmean(25,:) = mean([WPCA1_eva.models{i}.T(28,:);WPCA1_eva.models{i}.Tpred(7,:)]);
    WPCA1_eva.models{i}.Tmean(26,:) = mean([WPCA1_eva.models{i}.T(29,:);WPCA1_eva.models{i}.Tpred(8,:)]);
    WPCA1_eva.models{i}.Tmean(27:28,:) = WPCA1_eva.models{i}.T(30:31,:);
    WPCA1_eva.models{i}.Tmean(29,:) = mean([WPCA1_eva.models{i}.T(32,:);WPCA1_eva.models{i}.Tpred(9,:)]);
    WPCA1_eva.models{i}.Tmean(30:31,:) = WPCA1_eva.models{i}.T(33:34,:);
    WPCA1_eva.models{i}.Tmean(32,:) = mean([WPCA1_eva.models{i}.T(35,:);WPCA1_eva.models{i}.Tpred(10,:)]);
    WPCA1_eva.models{i}.Tmean(33,:) = WPCA1_eva.models{i}.T(36,:);
    
    % Computing variances of each PC
    
    X1 = var(WPCA1_eva.models{i}.T([1,2,15,26,37],:));
    X2 = var([WPCA1_eva.models{i}.T(3,:);WPCA1_eva.models{i}.Tpred(6,:)]);
    X3 = var([WPCA1_eva.models{i}.T(9,:);WPCA1_eva.models{i}.Tpred(1,:)]);
    X4 = var([WPCA1_eva.models{i}.T(13,:);WPCA1_eva.models{i}.Tpred(2,:)]);
    X5 = var([WPCA1_eva.models{i}.T(18,:);WPCA1_eva.models{i}.Tpred(3,:)]);
    X6 = var([WPCA1_eva.models{i}.T(20,:);WPCA1_eva.models{i}.Tpred(4,:)]);
    X7 = var([WPCA1_eva.models{i}.T(23,:);WPCA1_eva.models{i}.Tpred(5,:)]);
    X8 = var([WPCA1_eva.models{i}.T(28,:);WPCA1_eva.models{i}.Tpred(7,:)]);
    X9 = var([WPCA1_eva.models{i}.T(29,:);WPCA1_eva.models{i}.Tpred(8,:)]);
    X10 = var([WPCA1_eva.models{i}.T(32,:);WPCA1_eva.models{i}.Tpred(9,:)]);
    X11 = var([WPCA1_eva.models{i}.T(35,:);WPCA1_eva.models{i}.Tpred(10,:)]);
    
    for j = 1:size(X1,2);
        
        WPCA1_eva.models{i}.stdpooled(j) = sqrt(((4*X1(1,j))+X2(1,j)+X3(1,j)+X4(1,j)+X5(1,j)+...
            X6(1,j)+X7(1,j)++X8(1,j)+X9(1,j)+X10(1,j)+X11(1,j))/14);
        
    end
    
    t = tinv(0.95,14); % t-value at 95% conf. level
    
    for k = 1:size(X1,2);
    
        for j = 1:size(WPCA1_eva.models{i}.Tmean,1);
            WPCA1_eva.models{i}.Tconfint(j,k) = abs((WPCA1_eva.models{i}.Tmean(j,k)+(WPCA1_eva.models{i}.stdpooled(k)*t/sqrt(25))) -...
                (WPCA1_eva.models{i}.Tmean(j,k)-(WPCA1_eva.models{i}.stdpooled(k)*t/sqrt(25)))); % Computing confidence intervals
        end
    end 
    
    clear X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 i j k t
    
end

%% Plotting mean Scores and Confidence Intervals

i = 1; % selected model
pcx = 1;
pcy = 2;

samples = Data.info([1,3:9,11:14,16,18:20,22:23,25:27,29:30,33:34,36,38:40,42:44,46],1:2);

s1cal = cellfun('isempty', strfind(samples(:,2),'S'));
s2cal = find(s1cal == 0);
figure;errorbar(WPCA1_eva.models{i}.Tmean(s2cal,pcx),WPCA1_eva.models{i}.Tmean(s2cal,pcy),WPCA1_eva.models{i}.Tconfint(s2cal,pcy),'s','MarkerFaceColor','blue','Color','blue')
hold on
herrorbar(WPCA1_eva.models{i}.Tmean(s2cal,pcx),WPCA1_eva.models{i}.Tmean(s2cal,pcy),WPCA1_eva.models{i}.Tconfint(s2cal,pcx),'blue')
text(WPCA1_eva.models{i}.Tmean(s2cal,pcx),WPCA1_eva.models{i}.Tmean(s2cal,pcy),samplesNew(s2cal,1))

d1cal = cellfun('isempty', strfind(samples(:,2),'D'));
d2cal = find(d1cal == 0);
errorbar(WPCA1_eva.models{i}.Tmean(d2cal,pcx),WPCA1_eva.models{i}.Tmean(d2cal,pcy),WPCA1_eva.models{i}.Tconfint(d2cal,pcy),'s','MarkerFaceColor','green','Color','green')
herrorbar(WPCA1_eva.models{i}.Tmean(d2cal,pcx),WPCA1_eva.models{i}.Tmean(d2cal,pcy),WPCA1_eva.models{i}.Tconfint(d2cal,pcx),'green')
text(WPCA1_eva.models{i}.Tmean(d2cal,pcx),WPCA1_eva.models{i}.Tmean(d2cal,pcy),samplesNew(d2cal,1))


%legend('QC','S','D')

xlabel(['PC1 ',num2str(pcx),'(',num2str(round(WPCA1_eva.models{i}.SSQ(pcx,3),2)),'%',')'])
ylabel(['PC',num2str(pcy),'(',num2str(round(WPCA1_eva.models{i}.SSQ(pcy,3),2)),'%',')'])
title(['rt1D (min) > ',num2str(WPCA1_eva.rt1int(i))])

clear s1cal s2cal d1cal d2cal samples pcx 
%Creating the indexes in the unfolded data according to the evolutiopcy

%% PLOTTING 2D LOADINGS in the EVOLUTIONARY WPCA MODELS
sicpos = [9]; %position of the target sic
k = 1; %model to analyse
pc = 2;

a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(1)}(1,1):rtLim{WPCA1_eva.sics(1)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');
b = size(Xal{1,WPCA1_eva.sics(1)},1)*(size(Xal{1,WPCA1_eva.sics(1)},2)-a+1); 
ind2 = b; clear a b

for i = 2:length(WPCA1_eva.sics);
   
 a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(i)}(1,1):rtLim{WPCA1_eva.sics(i)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');

 if isempty(a) == 0;
    b = size(Xal{1,WPCA1_eva.sics(i)},1)*(size(Xal{1,WPCA1_eva.sics(i)},2)-a+1);
    ind2 = [ind2 ind2(i-1)+b];
 else
    ind2 = [ind2 ind2(i-1)];
 end
 
end

ind2 = [0 ind2]; clear a b i;

X2 = WPCA1_eva.models{k}.P(ind2(sicpos)+1:ind2(sicpos+1),pc);

a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(sicpos)}(1,1):rtLim{WPCA1_eva.sics(sicpos)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');

X2 = reshape(X2,size(Xal{1,WPCA1_eva.sics(sicpos)},1),(size(Xal{1,WPCA1_eva.sics(sicpos)},2)-a+1));

figure;imagesc(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(sicpos)}(1,1)+a:rtLim{WPCA1_eva.sics(sicpos)}(1,2),1),Batch2.rt2{1,2}(rtLim{WPCA1_eva.sics(sicpos)}(2,1):rtLim{WPCA1_eva.sics(sicpos)}(2,2)),X2)
axis xy
xlabel('retention time 1D (min)');
ylabel('retention time 2D (s)');

clear ind2 sicpos X2 pc a k      

%% PLOTTING THE LOADINGS FROM ALL THE EICs SIMULTANEOUSLY

A = zeros(rtLim{1,1}(2,2),rtLim{1,1}(1,2));
k = 1; %model to analyse
pc = 1;

for sicpos = 1:3;
%Creating the indexes in the unfolded data according to the evolutionary
%model k

a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(1)}(1,1):rtLim{WPCA1_eva.sics(1)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');
b = size(Xal{1,WPCA1_eva.sics(1)},1)*(size(Xal{1,WPCA1_eva.sics(1)},2)-a+1); 
ind2 = b; clear a b

for i = 2:length(WPCA1_eva.sics);
   
 a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(i)}(1,1):rtLim{WPCA1_eva.sics(i)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');

 if isempty(a) == 0;
    b = size(Xal{1,WPCA1_eva.sics(i)},1)*(size(Xal{1,WPCA1_eva.sics(i)},2)-a+1);
    ind2 = [ind2 ind2(i-1)+b];
 else
    ind2 = [ind2 ind2(i-1)];
 end
 
end

ind2 = [0 ind2]; clear a b i;

X2 = WPCA1_eva.models{k}.P(ind2(sicpos)+1:ind2(sicpos+1),pc);

a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(sicpos)}(1,1):rtLim{WPCA1_eva.sics(sicpos)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');

X2 = reshape(X2,size(Xal{1,WPCA1_eva.sics(sicpos)},1),(size(Xal{1,WPCA1_eva.sics(sicpos)},2)-a+1));

rt1int = Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(sicpos)}(1,1)-1+a:rtLim{WPCA1_eva.sics(sicpos)}(1,2),1);
rt2int = Batch2.rt2{1,2}(rtLim{WPCA1_eva.sics(sicpos)}(2,1):rtLim{WPCA1_eva.sics(sicpos)}(2,2));

A(find(Batch2.rt2{1,2} == rt2int(1,1)):find(Batch2.rt2{1,2} == rt2int(end,1)),find(Batch2.rt1{1,2} == rt1int(1,1)):find(Batch2.rt1{1,2} == rt1int(end,1))) = X2;
clear X2 rt1int rt2int

end

figure;imagesc(Batch2.rt1{1,2}(:,1),Batch2.rt2{1,2}(:,1),A)
axis xy
clear A pc k sicpos ind2
%% CREATING BAR PLOTS THROUGH LOAGINGS PEAKS OF PAHs

% Creating the indexes

sicpos = [22]; %position of the target sics
k = 7; %model to analyse
pc = 1;

a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(1)}(1,1):rtLim{WPCA1_eva.sics(1)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');
b = size(Xal{1,WPCA1_eva.sics(1)},1)*(size(Xal{1,WPCA1_eva.sics(1)},2)-a+1); 
ind2 = b; clear a b

for i = 2:length(WPCA1_eva.sics);
   
 a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(i)}(1,1):rtLim{WPCA1_eva.sics(i)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');

 if isempty(a) == 0;
    b = size(Xal{1,WPCA1_eva.sics(i)},1)*(size(Xal{1,WPCA1_eva.sics(i)},2)-a+1);
    ind2 = [ind2 ind2(i-1)+b];
 else
    ind2 = [ind2 ind2(i-1)];
 end
 
end

ind2 = [0 ind2]; clear b i;

b = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(sicpos)}(1,1):rtLim{WPCA1_eva.sics(sicpos)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');

X2 = WPCA1_eva.models{k}.P(ind2(sicpos)+1:ind2(sicpos+1),pc);
X3 = reshape(X2,size(Xal{1,WPCA1_eva.sics(sicpos)},1),(size(Xal{1,WPCA1_eva.sics(sicpos)},2)-b+1));


figure;plot(X2),hold on,plot(abs(X2)), clear ind2 b,

figure;imagesc(X3)
axis xy

c = find(X2 < 0);

WPCA1_eva.models{k}.PAHint{pc,sicpos} = [1,length(X2)];
a = WPCA1_eva.models{k}.PAHint{1,sicpos};
lim = 1E-3;
lim2 = 40;
figure; findpeaks(abs(X2(a(1):a(end),1)),'minpeakheight',lim, 'minpeakdistance',lim2);

%%
[~,l,~,p] = findpeaks(abs(X2(a(1):a(end),1)),'minpeakheight',lim, 'minpeakdistance',lim2);

clear X3 c a

c = find(X2(l,1) < 0); 

if isempty(c) == 1;
    WPCA1_eva.models{k}.PAH(pc,sicpos) = sum(p);
else
    d = p(c).*-1;
    
    WPCA1_eva.models{k}.PAH(pc,sicpos) = sum(d);
end

clear d pc k l p X2 sicpos lim lim2

%% Plotting the bars
k = 1;
figure; bar(10:30,WPCA1_eva.models{k}.PAH(1,10:end)), axis tight
hold on
bar(10:30,WPCA1_eva.models{k}.PAH(2,10:end),'red'), axis tight


%% WITHIN NORMALIZATION  of SPECIFIC SICs

%% NAPHTHENES

WPCANap.sic = 1:5;

k = 5; %(rt1D > 37.18 min)

a = find(Batch2.rt1{1,2}(rtLim{WPCANap.sic(1)}(1,1):rtLim{WPCANap.sic(1)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');
b = size(Xal{1,WPCANap.sic(1)},1)*(size(Xal{1,WPCANap.sic(1)},2)-a+1); 
ind2 = b; clear a b

for i = 2:length(WPCANap.sic);
   
 a = find(Batch2.rt1{1,2}(rtLim{WPCANap.sic(i)}(1,1):rtLim{WPCANap.sic(i)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');

 if isempty(a) == 0;
    b = size(Xal{1,WPCANap.sic(i)},1)*(size(Xal{1,WPCANap.sic(i)},2)-a+1);
    ind2 = [ind2 ind2(i-1)+b];
 else
    ind2 = [ind2 ind2(i-1)];
 end
 
end

ind2 = [0 ind2]; clear b i;

WPCANap.ind2 = ind2; clear ind2 a
%%
Xmodel = [];
sicpos = 1:5;
k = 5;
for j = 1:length(WPCANap.sic);
    
    a = find(Batch2.rt1{1,2}(rtLim{WPCANap.sic(j)}(1,1):rtLim{WPCANap.sic(j)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');

    for l = 1:size(Xal3,1);
        A = Xal3{l,sicpos(j)}(:,a:end);
        X(l,:) = reshape(A,1,size(A,1)*size(A,2));
    end
    %X = normr(X);  % normalizing each SIC individually
    
    Xmodel = [Xmodel X];
    
    clear X 
end  

WPCANap.Data = normr(Xmodel);
clear Xmodel sicpos

%% COMPUTING WPCA MODEL

%Extracting the corresponding weights
   
X1 = WPCANap.Data([2,1,17,32,47],:);
    
w = 1./(std(X1)./mean(X1));

c = find(isnan(w) == 1);
w(1,c) = zeros(1,length(c));

clear X1 c
        
        
% Calculating the models
    
for k = 1:size(WPCANap.Data,1);
    X(k,:) = WPCANap.Data(k,:).*w; %downscaling variables
end
    
Xcal = X(WPCA1_eva.cal,:);
Xval = X(WPCA1_eva.val,:);
    
[~,~,rmsecv,~,~,~,~] = crossval(Xcal,[],'pca',{'vet' 5},5,optcv);
a = find(rmsecv == min(rmsecv),1,'first');
    
model = pca(Xcal,a,optpca);
    
WPCANap.T = model.loads{1,1}; % extracting scores
    
for k = 1:size(model.loads{2,1},2);
    WPCANap.P(:,k) = model.loads{2,1}(:,k)./w';
    c = find(isnan(WPCANap.P(:,k)) == 1);
    WPCANap.P(c,k) = zeros(length(c),1);
    c = find(isfinite(WPCANap.P(:,k)) == 0);
    WPCANap.P(c,k) = zeros(length(c),1);
end
  
WPCANap.SSQ = model.detail.ssq(1:a,:);
    
pred = pca(Xval,model,optpca);
WPCANap.Tpred = pred.loads{1,1}; % extracting validation scores
    
clear model pred k Xcal Xval a rmsecv w X c

%% Plotting the PC1 x PC2 scores (with validation)
i = 5;

qc = cellfun('isempty', strfind(Data.info(WPCA1_eva.cal,2),'QC'));
qc2 = find(qc == 0);
A = Data.info(WPCA1_eva.cal,1);
figure;plot(WPCANap.T(qc2,1),WPCANap.T(qc2,2),'o','MarkerFaceColor','red','Color','red')
text(WPCANap.T(qc2,1),WPCANap.T(qc2,2),strtok(A(qc2,1),'.'))

hold on

s1cal = cellfun('isempty', strfind(Data.info(WPCA1_eva.cal,2),'S'));
s2cal = find(s1cal == 0);
plot(WPCANap.T(s2cal,1),WPCANap.T(s2cal,2),'s','MarkerFaceColor','blue','Color','blue')
text(WPCANap.T(s2cal,1),WPCANap.T(s2cal,2),strtok(A(s2cal,1),'.'))

d1cal = cellfun('isempty', strfind(Data.info(WPCA1_eva.cal,2),'D'));
d2cal = find(d1cal == 0);
plot(WPCANap.T(d2cal,1),WPCANap.T(d2cal,2),'s','MarkerFaceColor','green','Color','green')
text(WPCANap.T(d2cal,1),WPCANap.T(d2cal,2),strtok(A(d2cal,1),'.'))

B = Data.info(WPCA1_eva.val,1);
s1val = cellfun('isempty', strfind(Data.info(WPCA1_eva.val,2),'S'));
s2val = find(s1val == 0);
plot(WPCANap.Tpred(s2val,1),WPCANap.Tpred(s2val,2),'*','MarkerFaceColor','blue','Color','blue')
text(WPCANap.Tpred(s2val,1),WPCANap.Tpred(s2val,2),strtok(B(s2val,1),'.'))

d1val = cellfun('isempty', strfind(Data.info(WPCA1_eva.val,2),'D'));
d2val = find(d1val == 0);
plot(WPCANap.Tpred(d2val,1),WPCANap.Tpred(d2val,2),'*','MarkerFaceColor','green','Color','green')
text(WPCANap.Tpred(d2val,1),WPCANap.Tpred(d2val,2),strtok(B(d2val,1),'.'))

legend('QC','S','D','Sval','Dval')

xlabel(['PC1 ','(',num2str(round(WPCANap.SSQ(1,3),2)),'%',')'])
ylabel(['PC2 ','(',num2str(round(WPCANap.SSQ(2,3),2)),'%',')'])
title(['rt1D (min) > ',num2str(WPCA1_eva.rt1int(i))])
clear qc qc2 s1cal s2cal d1cal d2cal s1val s2val d1val d2val A B i


%% PLOTTING 2D LOADINGS in the EVOLUTIONARY WPCA MODELS
sicpos = [4]; %position of the target sics
k = 1; %model to analyse
pc = 2;


%Creating the indexes in the unfolded data according to the evolutionary
%model k

a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(1)}(1,1):rtLim{WPCA1_eva.sics(1)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');
b = size(Xal{1,WPCA1_eva.sics(1)},1)*(size(Xal{1,WPCA1_eva.sics(1)},2)-a+1); 
ind2 = b; clear a b

for i = 2:length(WPCA1_eva.sics);
   
 a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(i)}(1,1):rtLim{WPCA1_eva.sics(i)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');

 if isempty(a) == 0;
    b = size(Xal{1,WPCA1_eva.sics(i)},1)*(size(Xal{1,WPCA1_eva.sics(i)},2)-a+1);
    ind2 = [ind2 ind2(i-1)+b];
 else
    ind2 = [ind2 ind2(i-1)];
 end
 
end

ind2 = [0 ind2]; clear a b i;

%X2 = WPCA1_eva.models{k}.P(ind2(sicpos)+1:ind2(sicpos+1),pc);

X2 = WPCANap.P(ind2(sicpos)+1:ind2(sicpos+1),pc);
a = find(Batch2.rt1{1,2}(rtLim{WPCA1_eva.sics(sicpos)}(1,1):rtLim{WPCA1_eva.sics(sicpos)}(1,2),1) >= WPCA1_eva.rt1int(k),1,'first');

X2 = reshape(X2,size(Xal{1,WPCA1_eva.sics(sicpos)},1),(size(Xal{1,WPCA1_eva.sics(sicpos)},2)-a+1));

figure; imagesc(X2);axis xy

%figure;imagesc(Batch2.rt1{1,2}(rtLim{WPCA2.sics(sicpos)}(1,1)+a:rtLim{WPCA2.sics(sicpos)}(1,2),1),Batch2.rt2{1,2}(rtLim{WPCA2.sics(sicpos)}(2,1):rtLim{WPCA2.sics(sicpos)}(2,2)),X2)
%axis xy

clear ind2 sicpos X2 pc a k      




%% Assigning DIESELS and SPILLS in WPCA2_eva models (KNOWN cases)

WPCA1_eva.Cases1.c1.samples = {'3183-1';'3183-2'};
WPCA1_eva.Cases1.c2.samples = {'3221-1';'3221-2'};
WPCA1_eva.Cases1.c3.samples = {'J3253-8';'J3253-12';'J3253-13';'J3253-19';'J3253-20';'J3253-23'};
WPCA1_eva.Cases1.c4.samples = {'3243-4';'3243-7';'3243-8';'3243-11';'3243-12'};

samples = Data.info([1,3:9,11:14,16,18:20,22:23,25:27,29:30,33:34,36,38:40,42:44,46],1:2);

% Calculating Probabilities

for l = 1:length(WPCA1_eva.models);
    
    %WPCA1_eva.models{l}.Cases1.c4.Tmean = WPCA1_eva.models{l}.Tmean([30,29,28,27,26,25],:);
    
    for i = 1:length(WPCA1_eva.models{l}.stdpooled); % PCs
    
        k = 1; %spill sample in WPCA1_eva.models.Cases1.c4.samples
    
       
        for j = 2:5; %sources
            %q{j}(k,i) = abs(WPCA1_eva.models{l}.Cases1.c4.Tmean(k,i) - WPCA1_eva.models{l}.Cases1.c4.Tmean(j,i));
            q{j}(k,i) = abs(WPCA1_eva.models{l}.Cases1.c4.Tmean(k,i) - WPCA1_eva.models{l}.Cases1.c4.Tmean(j,i))/(WPCA1_eva.models{l}.stdpooled(i)*sqrt(2));
            WPCA1_eva.models{l}.Cases1.c4.prob{k}(i,j) = 2*(1 - tcdf(q{j}(k,i),14));
        
        end
    
    end
    for j = 2:5;
        %WPCA1_eva.models{l}.Cases1.c4.Tdist(k,j) = sum(q{j}(k,1:length(WPCA1_eva.models{l}.stdpooled))'.*(WPCA1_eva.models{l}.SSQ(1:length(WPCA1_eva.models{l}.stdpooled),3)./100))/(WPCA1_eva.models{l}.SSQ(length(WPCA1_eva.models{l}.stdpooled),4)/100);
        WPCA1_eva.models{l}.Cases1.c4.meanprob(k,j) = sum(WPCA1_eva.models{l}.Cases1.c4.prob{k}(:,j).*(WPCA1_eva.models{l}.SSQ(1:length(WPCA1_eva.models{l}.stdpooled),3)./100))/(WPCA1_eva.models{l}.SSQ(length(WPCA1_eva.models{l}.stdpooled),4)/100);

    end
    clear i j k
end

clear q l
%% Plotting the Distances between Spills and Sources
hold on
k = 6; %source

for i = 1:length(WPCA1_eva.models);
    for j = 1:3; %spill
        a(j,i) = mean(WPCA1_eva.models{i}.Cases1.c3.Tdist(j,k),2);    
    end
    b(i) = sum(WPCA1_eva.models{i}.Tconfint(1,:)'.*(WPCA1_eva.models{i}.SSQ(1:length(WPCA1_eva.models{i}.stdpooled),3)./100))/(WPCA1_eva.models{i}.SSQ(length(WPCA1_eva.models{i}.stdpooled),4)/100);
end

%figure;
for k = 1:3; %spill
    plot([1:length(WPCA1_eva.models)],a(k,:),'s','Color','blue','MarkerFaceColor','blue','LineStyle','-'), hold on,
    errorbar([1:length(WPCA1_eva.models)],a(k,:),b,'Color','blue')
end

clear a b i j k

%%

for i = 1:length(WPCA1_eva.models);
    for k = 2:5;
        a(k,i) = mean(WPCA1_eva.models{i}.Cases1.c4.Tdist(1,k),2);    
    end
    b(i) = sum(WPCA1_eva.models{i}.Tconfint(1,:)'.*(WPCA1_eva.models{i}.SSQ(1:length(WPCA1_eva.models{i}.stdpooled),3)./100))/(WPCA1_eva.models{i}.SSQ(length(WPCA1_eva.models{i}.stdpooled),4)/100);
end


figure;
for t = 2:size(a,1);
    plot([1:length(WPCA1_eva.models)],a(t,:),'s','Color','green','MarkerFaceColor','blue','LineStyle','-'), hold on,
    errorbar([1:length(WPCA1_eva.models)],a(t,:),b,'Color','blue')
end
clear a b i j k t



%% Assigning DIESELS and SPILLS in WPCA2_eva models (Unknown cases)
samples2 = Data.info([1,3:9,11:14,16,18:20,22:23,25:27,29:30,33:34,36,38:40,42:44,46],1:2);

samples = samples2([2:9,11:15,17:21,23:27,29:33],:);

spill = cellfun('isempty', strfind(samples(:,2),'S'));
spill2 = find(spill == 0);
diesel = cellfun('isempty', strfind(samples(:,2),'D'));
diesel2 = find(diesel == 0);

  
 WPCA1_eva.Cases2.samplesSpill = samples(spill2,:);
 WPCA1_eva.Cases2.samplesDiesel = samples(diesel2,:);
    

clear spill spill2 diesel diesel2 samples samples2


%% Calculating Scores distances

samples2 = Data.info([1,3:9,11:14,16,18:20,22:23,25:27,29:30,33:34,36,38:40,42:44,46],1:2);
samples = samples2([2:9,11:15,17:21,23:27,29:33],:);

spill = find(cellfun('isempty', strfind(samples(:,2),'S')) == 0);
diesel = find(cellfun('isempty', strfind(samples(:,2),'D')) == 0);


for k  = 1:length(WPCA1_eva.models);

    T = WPCA1_eva.models{k}.Tmean([2:9,11:15,17:21,23:27,29:33],:);
    
    for i = 1:length(WPCA1_eva.Cases2.samplesSpill);
    
        for j = 1:length(WPCA1_eva.Cases2.samplesDiesel);
        
            Tmean = T([spill(i),diesel(j)],:);
        
            for l = 1:size(Tmean,2);
                q(i,j,l) = abs(Tmean(1,l) - Tmean(2,l));
                %q = (abs(Tmean(1,l) - Tmean(2,l)))/(WPCA1_eva.models{k}.stdpooled(l)*sqrt(2));
                %WPCA1_eva.models{k}.Cases2.prob(i,j,l) = 2*(1 - tcdf(q,14));
            end
            
            
            a = reshape(q(i,j,:),1,size(q(i,j,:),3));
            %a = reshape(WPCA1_eva.models{k}.Cases2.prob(i,j,:),1,size(WPCA1_eva.models{k}.Cases2.prob(i,j,:),3));
            WPCA1_eva.models{k}.Cases2.Tdist(i,j) = sum(a'.*(WPCA1_eva.models{k}.SSQ(1:size(T,2),3)./100))/(WPCA1_eva.models{k}.SSQ(size(T,2),4)/100);
            %WPCA1_eva.models{k}.Cases2.meanprob(i,j) = sum(a'.*(WPCA1_eva.models{k}.SSQ(1:size(T,2),3)./100))/(WPCA1_eva.models{k}.SSQ(size(T,2),4)/100);

            clear a l
       end
    clear j
           
    end

    clear i q
end
clear diesel spill samples samples2 k T

    
   

%% Calculating Probability of belonging to one source

samples2 = Data.info([1,3:9,11:14,16,18:20,22:23,25:27,29:30,33:34,36,38:40,42:44,46],1:2);
samples = samples2([2:9,11:15,17:21,23:27,29:33],:);

spill = find(cellfun('isempty', strfind(samples(:,2),'S')) == 0);
diesel = find(cellfun('isempty', strfind(samples(:,2),'D')) == 0);


for k  = 1:length(WPCA1_eva.models);

    T = WPCA1_eva.models{k}.Tmean([2:9,11:15,17:21,23:27,29:33],:);
    
    for i = 1:length(WPCA1_eva.Cases2.samplesSpill);
    
        for j = 1:length(WPCA1_eva.Cases2.samplesDiesel);
        
            Tmean = T([spill(i),diesel(j)],:);
        
            for l = 1:size(Tmean,2);
                %q(i,j,l) = abs(Tmean(1,l) - Tmean(2,l));
                q = (abs(Tmean(1,l) - Tmean(2,l)))/(WPCA1_eva.models{k}.stdpooled(l)*sqrt(2));
                WPCA1_eva.models{k}.Cases2.prob(i,j,l) = 2*(1 - tcdf(q,14));
            end
            
            
            %a = reshape(q(i,j,:),1,size(q(i,j,:),3));
            a = reshape(WPCA1_eva.models{k}.Cases2.prob(i,j,:),1,size(WPCA1_eva.models{k}.Cases2.prob(i,j,:),3));
            %WPCA1_eva.models{k}.Cases2.Tdist(i,j) = sum(a'.*(WPCA1_eva.models{k}.SSQ(1:size(T,2),3)./100))/(WPCA1_eva.models{k}.SSQ(size(T,2),4)/100);
            WPCA1_eva.models{k}.Cases2.meanprob(i,j) = sum(a'.*(WPCA1_eva.models{k}.SSQ(1:size(T,2),3)./100))/(WPCA1_eva.models{k}.SSQ(size(T,2),4)/100);

            clear a l
       end
    clear j
           
    end

    clear i q
end
clear diesel spill samples samples2 k T

    


%% PERFORMING OIL-SOURCE CORRELATION AFTER ELIMINATING HIGHY WEATHERED SPILLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computing new W now considering high degree of weathering (evaporation)
% Samples 3242-3a 3242-3b 11781 3210-1 (high evaporation)



















%%
X(:,:,1) = sum(Batch3.SIC_GPS{2},3);
X(:,:,2) = Batch3.SIC_GPS{1}(:,:,i);
X(:,:,3) = Batch3.SIC_GPS{2}(:,:,i);
X(:,:,4) = Batch3.SIC_GPS{1}(:,:,i);
X(:,:,5) = Batch3.SIC_GPS{1}(:,:,i);

%%
figure;imagesc(Batch3.rt1{13},Batch3.rt2{13},sum(Batch3.SIC_GPS{13},3))
axis xy

%%
i = 32;
    
figure;imagesc(Batch1.SIC_GPSsub{3,i})
figure;imagesc(Batch2.SIC_GPSsub{3,i})
figure;imagesc(Batch3.SIC_GPSsub{12,i})
figure;imagesc(Batch2.SIC_GPSsub{2,i})
%figure;imagesc(Batch1.rt1{i},Batch1.rt2{i},sum(retX(:,:,i),3))

%axis xy

%%


imagesc(sum(Batch3.SIC_GPS{13}(:,:,1),3),size(Batch3.SIC_GPS{13},1)*size(Batch3.SIC_GPS{13},2),1)

%end
%% Handling Figures

ah = gca;
get(2,'Children');
ah(2) = ans(1);
ah(1).Parent;
linkaxes(ah,'xy');
d = sum((x-y).^2).^0.5;
WPCA2.Cases2.c1.Tmean = WPCA2.Tmean([12,23,32,15],:);
WPCA2.Cases2.c2.Tmean = WPCA2.Tmean([6,2,14],:);
WPCA2.Cases2.c3.Tmean = WPCA2.Tmean([7,2,14],:);
WPCA2.Cases2.c4.Tmean = WPCA2.Tmean([8,2,14],:);

WPCA1_eva.models{k}.Cases2.prob(i,j,l) = 2*(1 - tcdf(q,14));
%% Saving the data
save('Data_RCbck.mat');
%%
find(out.mass_values >= 212.05 & out.mass_values < 212.1);
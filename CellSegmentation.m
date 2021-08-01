%%%%% Cell segmentation
%%%%% Sanjeev Uthishtran - July 2021
%%%%% Updated - 22/07/2021

close all
clear integratedGrayValues

[fileA,pathA] = uigetfile('*.lsm','SELECT MASTER IMAGE','SELECT MASTER IMAGE');
disp('Master Image Selected')
DESTINATION_FOLDER = uigetdir('C:\','SELECT DESTINATION FOLDER');
disp('Destination folder selected')
%%%%Cell Channel Stitching
[Data] = lsmread(fullfile(pathA,fileA));
holi = size(Data);
totalpixels = holi(2)*holi(3)*holi(4)*holi(5);
channel2 = Data([2:2:totalpixels]);
C2mod = reshape(channel2,[holi(4),holi(5),holi(3)]);

StitchedvolumeA = [];
PolA = [];
t = 1;
for t = 1:holi(3)
    Slicepixels = C2mod([t:holi(3):totalpixels/2]);
    FinalsliceA = reshape(Slicepixels,[holi(4),holi(5)]);
    PolA = [PolA,sum(FinalsliceA,'all')];
    StitchedvolumeA = cat(3,StitchedvolumeA,FinalsliceA);
end
MidA = find(PolA == max(PolA));
A = sum(StitchedvolumeA(:,:,MidA-5:MidA+5),3);

%%%%Nuclei Channel Stitching
channel1 = Data([1:2:totalpixels]);
C1mod = reshape(channel1,[holi(4),holi(5),holi(3)]);

StitchedvolumeB = [];
PolB = [];
t = 1;
for t = 1:holi(3)
    SlicepixelsB = C1mod([t:holi(3):totalpixels/2]);
    FinalsliceB = reshape(SlicepixelsB,[holi(4),holi(5)]);
    PolB = [PolB,sum(FinalsliceB,'all')];
    StitchedvolumeB = cat(3,StitchedvolumeB,FinalsliceB);
end
MidB = find(PolB == max(PolB));
B = sum(StitchedvolumeB(:,:,MidB-3:MidB+3),3);
%%%%%%%%%%
disp('Generating original master image')
figure, imshow(A), colorbar,caxis([0 2000]),title('Original Master Image')
saveas(figure(1),fullfile(DESTINATION_FOLDER,'1_Original_Master_Image.png'))
%%%%%%%%%%
disp('Generating original master image')
figure, imshow(B), colorbar,caxis([0 2000]),title('Original Nuclei Image')
saveas(figure(2),fullfile(DESTINATION_FOLDER,'2_Nuclei_Master_Image.png'))

%Option for filters
prompt = 'Do you want to apply filters? Y/N: ';
x = input(prompt,'s');

%Condition for filters or not
if x == 'Y'
    prompt = {'Area Threshold:','Threshold Adjusting Factor:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'4000','100'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    
    Areathresh = str2num(answer{1});
    Adjustfactor = str2num(answer{2});

    MedA = medfilt2(A,[4 4]);
    
    MaxA = nlfilter(MedA,[4 4],@(x)max(x(:)));
    
    Thresh = MaxA;
    CustomMed = median(MaxA,'all')+Adjustfactor;
    Thresh(Thresh<CustomMed) = 0;
    Thresh(Thresh>=CustomMed) = 1;

    flams = imfill(Thresh,'holes');
    seDA = strel('diamond',1);
    BWfinalB = imerode(flams, seDA);
    BWfinalB = imerode(BWfinalB, seDA);

elseif x == 'N'
    prompt = {'Area Threshold:','Edge Sensitivity:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'4000','0.15'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    Areathresh = str2num(answer{1});
    [~,thresholdB] = edge(A, 'sobel');
    
    fudgeFactorB = str2num(answer{2});
    
    Thresh = edge(A,'sobel',thresholdB*fudgeFactorB);
    
    se90A = strel('line',3,90);
    se0A = strel('line',3,0);
    BWsdilB = imdilate(Thresh,[se90A,se0A]);
    BWdfillB = imfill(BWsdilB, 'holes');
    
    seDA = strel('diamond',1);
    BWfinalB = imerode(BWdfillB, seDA);
    BWfinalB = imerode(BWfinalB, seDA);
else
    disp('Invalid input: Please run the code again')
    return
end
%%%%%%%%%%

%%%%%%%%%%
disp('Generating segmented master image')
figure, imshow(BWfinalB), colorbar,caxis([0 1]),title('Segmented Master Image')
saveas(figure(3),fullfile(DESTINATION_FOLDER,'3_Segmented_Master_Image.png'))


%Option for connected sample
prompt = 'Are your cells connected after segmentation? Y/N: ';
y = input(prompt,'s');

if y == 'Y'
    prompt = {'Segmentation Factor:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'2'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    
    Segmentationfactor = str2num(answer{1});

   
    Sobelsegmentation = BWfinalB;
    segmentopen = ~bwareaopen(~Sobelsegmentation,10);
    %%%
    D = -bwdist(~Sobelsegmentation);
    %%%
    disp('First Watershed')
    FirstWS = watershed(D);
    %%%
    segmentopen = Sobelsegmentation;
    segmentopen(FirstWS == 0) = 0;
    %%%
    mask = imextendedmin(D,Segmentationfactor);
    %%%
    disp('Second Watershed')
    D2 = imimposemin(D,mask);
    Ld2 = watershed(D2);
    bw3 = Sobelsegmentation;
    bw3(Ld2 == 0) = 0;
    BWfinalB = bw3;
end 


%%%%%%%% Finding objects in the segment and filtering using area

CC = bwconncomp(BWfinalB);
L = labelmatrix(CC);
stats = regionprops('table',L,'Area');

numstats = table2array(stats);
label = find(numstats(:)>Areathresh);
f = transpose(1:length(label));
AMI = numstats(label,:);
combine = [f,AMI];

areas1 = ismember(L,label);

disp('Generating Filtered Segments')
figure,imshow(areas1),colorbar,caxis([0 1]),title('Filtered Segments')
saveas(figure(4),fullfile(DESTINATION_FOLDER,'4_Filtered_Segments.png'))


%%%%% Generate Integrated Intensities
[labeledImage, numCells] = bwlabel(areas1);
props = regionprops(labeledImage, A, 'PixelValues');
for k = 1 : numCells
    thisCellValues = props(k).PixelValues;
    integratedGrayValues(k) = sum(thisCellValues);
end

%%%%% Generating centroid stats for labelling diagram
finalstats = regionprops('table',labeledImage,A,'Centroid','Area','MeanIntensity');
finalstatsarr = table2array(finalstats);

x = finalstatsarr(:,2);
y = finalstatsarr(:,3);

disp('Generating Filtered Segments')
figure,imshow(labeledImage),title('Labeled segments')


hold on
for ii = 1:length(x)
    text(x(ii),y(ii),num2str(ii),'Color','m','FontSize',14)
end

saveas(figure(5),fullfile(DESTINATION_FOLDER,'5_Labeled_Segments.png'))

%%%%%%%%

disp('Generating labelled segmented region over brightfield')
figure
imshow(A),colorbar,caxis([0 2000])
hold on        
        for ii = 1:length(x)
             text(x(ii),y(ii),num2str(ii),'Color','m','FontSize',14)
        end
[Bla,Lla] = bwboundaries(areas1,'noholes');
for k = 1:length(Bla)
   boundary = Bla{k};
   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1)
end
saveas(figure(6),fullfile(DESTINATION_FOLDER,'6_Labelled_Segmented_Region_Over_Brightfield.png'))

%%%%%%%%

disp('Generating labelled segmented region over brightfield')
figure
imshow(A),colorbar,caxis([0 2000])
hold on        
        for ii = 1:length(x)
             text(x(ii),y(ii),num2str(ii),'Color','m','FontSize',14)
        end

for k = 1:length(Bla)
   boundary = Bla{k};
   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1)
end
saveas(figure(7),fullfile(DESTINATION_FOLDER,'7_Labelled_Segmented_Region_Over_Brightfield.png'))

%%%%%%%


%%%%%%%% Final table with stats
disp('Finalising statistics')
g = transpose(1:length(label));
FINALTABLE1 = table(g,finalstatsarr(:,1),finalstatsarr(:,4),transpose(integratedGrayValues),'VariableNames',{'Label','Area(Px)','MeanIntensity','Integrated Intensity'});
writetable(FINALTABLE1,fullfile(DESTINATION_FOLDER,'Statistics.csv'));

%%%%%%%%
figure
subplot(2,3,1),imshow(A), colorbar,caxis([0 5000]),title('Original Master Image')
subplot(2,3,2),imshow(B), colorbar,caxis([0 25000]),title('Original Nuclei Image')
subplot(2,3,3), imshow(BWfinalB), colorbar,caxis([0 1]),title('Segmented Master Image')
subplot(2,3,4),imshow(labeledImage),title('Filtered segments')
subplot(2,3,5),imshow(Segout),colorbar,caxis([0 50000]),title('Segmented Outline Over Original')
subplot(2,3,6),imshow(SegoutB),colorbar,caxis([0 50000]),title('Segmented Outline Over Nuclei Stain')
saveas(figure(8),fullfile(DESTINATION_FOLDER,'8_Compiled_Plots.png'))

%%%%%%%%

Puncta = [];

for ind = 1:length(label)
    chim = labeledImage;
    OGA = A;
    chim(labeledImage~=ind) = 0;
    loglabel = im2bw(chim);
    OGA(~loglabel) = 0;
    
    [~,thresholdB] = edge(OGA, 'sobel');
    fudgeFactorB = 6;
    BWsB = edge(OGA,'sobel',thresholdB*fudgeFactorB);

    se90A = strel('line',3,90);
    se0A = strel('line',3,0);

    BWsdilB = imdilate(BWsB,[se90A,se0A]);
    BWdfillB = imfill(BWsdilB, 'holes');

    seDA = strel('diamond',1);
    BWfinalB = imerode(BWdfillB, seDA);
    BWfinalB = imerode(BWfinalB, seDA);
    
    
    FF = bwconncomp(BWfinalB);
    P = labelmatrix(FF);
    numObjects = max(P(:));
    Puncta = [Puncta,numObjects];
    
end

FINALTABLE = table(g,finalstatsarr(:,1),finalstatsarr(:,4),transpose(integratedGrayValues),transpose(Puncta),'VariableNames',{'Label','Area(Px)','MeanIntensity','Integrated Intensity','Number of Puncta'});

disp('End')

    


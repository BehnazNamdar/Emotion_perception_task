clear all;
close all;

%settings
subjectNumber= input('subject code') ;
instructions= 'Press H for happy and S for sad faces';
fixationtime=.5;
stimulusduration= [.05 .1 .15 .2];
maskduration=.2;
delay=1;

escapeKey = KbName('ESCAPE');
sadKey = KbName('s');
happyKey = KbName('h');
nblock=10;
noise= [0 30 60 100];

%load images
path_ = cd;
path_ = [path_ '/images'];
filenames = dir(path_);
filenames = filenames(3:end);


%make noisy images
img= [];
image = [];

for a= 1 : size(filenames,1)
    
    image(a).name = num2str(a);
    image (a).img= imread(strcat(path_,'/',filenames(a).name));
    image(a).img = imresize (image(a).img , [300,300]);
    image(a).img = rgb2gray (image(a).img);
    imageid=a;  % 1 for sad , 2 for happy
    
    for b= 1:length(noise)
        
        k=length(noise);
        img((k*(a-1))+b).img= imnoise(image(a).img, 'salt & pepper', noise(b)/200);
        img((k*(a-1))+b).noiselevel= noise(b);
        img((k*(a-1))+b).id=imageid;
    end
    
end


% experimental matrix

for c=1:size(filenames,1)*length(noise)
   img(c).stimnumber = c;
end

stim_matrix = [img(:).stimnumber ;img(:).id; img(:).noiselevel]';
sd= vertcat(repmat(stimulusduration(1),size(filenames,1)*length(stimulusduration),1),...
    repmat(stimulusduration(2),size(filenames,1)*length(stimulusduration),1),...
    repmat(stimulusduration(3),size(filenames,1)*length(stimulusduration),1),...
    repmat(stimulusduration(4),size(filenames,1)*length(stimulusduration),1));
stim_matrix = horzcat(repmat(stim_matrix, length(stimulusduration), 1),sd);
exp_matrix = repmat (stim_matrix, nblock, 1);
exp_matrix = Shuffle(exp_matrix );
ntrials    = size(exp_matrix,1);

responsemat= zeros(ntrials , 3);

Screen('Preference', 'SkipSyncTests', 1);
[wPtr,rect]=Screen('Openwindow',max(Screen('Screens')),[128 128 128 ],[0 0 900 800]);

% fixation cross
centerx = rect(3)/2;
centery = rect(4)/2;
fixationhorz= [centerx-10 centery centerx+10 centery]';
fixationvert= [centerx centery-10 centerx centery+10]';

mask1 = imread(strcat(cd ,'/mask/m1.jpeg') );
mask1 = rgb2gray (imresize ( mask1 , [300,300]) );
mask2 = imread( strcat(cd ,'/mask/m2.jpeg') );
mask2 = rgb2gray (imresize ( mask2 , [300,300]) );

%experimental loop
for i = 1:ntrials
    
    if i==1
        
        DrawFormattedText(wPtr, instructions, 'center', 'center', 1);
        Screen('Flip', wPtr);
        KbStrokeWait;
        
    end
    
    Screen('DrawLine',wPtr,255,fixationhorz(1),fixationhorz(2),fixationhorz(3),fixationhorz(4),2);
    Screen('DrawLine',wPtr,255,fixationvert(1),fixationvert(2),fixationvert(3),fixationvert(4),2);
    Screen('Flip',wPtr);
    WaitSecs(fixationtime);
    
    mask1_= Screen('MakeTexture', wPtr ,mask1 ); 
    Screen('DrawTexture',wPtr,mask1_);
    Screen('Flip',wPtr);
    WaitSecs(maskduration);
    
    stimulustexture=Screen('MakeTexture',wPtr,img(exp_matrix(i,1)).img); 
    Screen ('DrawTexture',wPtr,stimulustexture);   
    Screen('Flip',wPtr);
    WaitSecs (exp_matrix(i,4));
    Screen('Flip',wPtr);
    
    mask2_= Screen('MakeTexture', wPtr ,mask2 ); 
    Screen('DrawTexture',wPtr,mask2_);
    Screen('Flip',wPtr);
    WaitSecs(maskduration);
    
    Screen('DrawLine',wPtr,255,fixationhorz(1),fixationhorz(2),fixationhorz(3),fixationhorz(4),2);
    Screen('DrawLine',wPtr,255,fixationvert(1),fixationvert(2),fixationvert(3),fixationvert(4),2);
    Screen('Flip',wPtr);
    WaitSecs(delay);
    
    respToBeMade = true;
    starttime = GetSecs;
    while respToBeMade
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(escapeKey) 
            sca;
            return
        elseif keyCode(sadKey)
            response = 1;
            respToBeMade = false;
        elseif keyCode(happyKey)
            response = 2;
            respToBeMade = false;
        elseif keyCode(KbName('space'))
            pause;
        end
    end
    endtime = GetSecs;
    rt = endtime - starttime;
    

    responsemat(i,1)=exp_matrix(i,2);  %truth
    responsemat(i,2)=exp_matrix(i,3);  %noise
    responsemat(i,3)=exp_matrix(i,4);
    responsemat(i,4)=response;
    responsemat(i,5)=rt;
    
end

sca;
save responsemat responsemat;

clear all;

%% Load data_set

path_ = cd;
path_ = [path_ '/data'];
filenames = dir(path_);
filenames = filenames(3:end);
ResponseMatrix =[];

for fi = 1 : size(filenames,1)
    load(strcat(path_ ,'/' , filenames(fi).name));
    ResponseMatrix(fi).name = fi;
end


for fi = 1 : size(filenames,1)
    load(strcat(path_ ,'/' , filenames(fi).name));
        ResponseMatrix(fi).true_rs = responsemat(: , 1);
        ResponseMatrix(fi).noise = responsemat(: , 2);
        ResponseMatrix(fi).stim_dur = responsemat(: , 3);
        ResponseMatrix(fi).rs = responsemat(: , 4);
        ResponseMatrix(fi).rt = responsemat(: , 5);
end



%% Performance and reaction time for each emotion per each display duration


perf_all_per_em_per_tim = [];
rt_all_per_em_per_tim   = [];


emo = categorical({'joy','sadness'});  % name of emotion to be displayed on x axis

emotions = unique([ResponseMatrix(1,:).true_rs]);
times      = unique([ResponseMatrix(1,:).stim_dur]);  

    
for si = 1 : size(ResponseMatrix,2)
 
    for ti = 1:length(times)
        
        for ci = 1:length(emotions)
            ind_h = [];
            ind_h = ismember([ResponseMatrix(1,si).true_rs],emotions(ci) ) & [ResponseMatrix(1,si).stim_dur] == times(ti) ;
            
            bhv = [];
            bhv = [ResponseMatrix(1,si).true_rs] == [ResponseMatrix(1,si).rs]  ;
            
            perf_all_per_em_per_tim(si,ci,ti) = 100*sum(bhv & ind_h)/sum(ind_h);
            rt_all_per_em_per_tim(si,ci,ti)   = nanmean(ResponseMatrix(1,si).rt(ind_h));
            
            
            
        end
    end
end


%% Plot psychometric function for each emotion on stimulus duration time

figure(1)
legends            = {'joy','sadness'};
for pm = 1:ci
    newcolors =[0    0.4470    0.7410
    0.8500    0.3250    0.0980];

    xlabel('Display duration (s)')
    ylabel('Mean Performance on all subjects (%)')

     plot(times, nanmean(squeeze(perf_all_per_em_per_tim(:,pm,:))),'LineWidth',4)
          hold on
      e=errorbar(times,nanmean(squeeze(perf_all_per_em_per_tim(:,pm,:))),nanstd(squeeze(perf_all_per_em_per_tim(:,pm,:)))/sqrt(size(perf_all_per_em_per_tim,1)),'HandleVisibility','off','Color',newcolors(pm,:),'LineWidth',1);
 
     legend(legends,'Location','northeast')
    xlim([0 0.25]);
    ylim([20 100]);
    title('Performance')
        set (gca, 'FontSize' ,15)
end
saveas( figure(1) , '1.jpg');

figure(2)
for pm=1:ci

    xlabel('Display diuration (s)')
    ylabel('Mean Reaction Time on all subjects (sec)')
     plot(times, nanmean(squeeze(rt_all_per_em_per_tim(:,pm,:))),'LineWidth',4)
        hold on

     errorbar(times,nanmean(squeeze(rt_all_per_em_per_tim(:,pm,:))),nanstd(squeeze(rt_all_per_em_per_tim(:,pm,:)))/sqrt(size(rt_all_per_em_per_tim,1)),'HandleVisibility','off','Color',newcolors(pm,:),'LineWidth',1)

    legend(legends,'Location','northeast')

    xlim([0 0.25]);
    title('Reaction Time')
    set (gca, 'FontSize' ,15)
end
saveas( figure(2) , '2.jpg');


%% Performance and reaction time for each emotion per each visual signal


perf_all_per_em_per_vs = [];
rt_all_per_em_per_vs   = [];
rt_all_per_em_per_vs_cr =[];
rt_all_per_em_per_vs_wr =[];
cl={};

emo = categorical({'joy','sadness'});  % name of emotion to be displayed on x axis

emotions = unique([ResponseMatrix(1,:).true_rs]);
times      = unique([ResponseMatrix(1,:).stim_dur]);  
noise = unique([ResponseMatrix(1,:).noise]);


for si = 1 : size(ResponseMatrix,2)
 
    for ni = 1:length(noise)
        
        for ci = 1:length(emotions)
            ind_h = [];
            ind_h = ismember([ResponseMatrix(1,si).true_rs],emotions(ci) ) & [ResponseMatrix(1,si).noise] == noise(ni) ;
            
            bhv = [];
            bhv = [ResponseMatrix(1,si).true_rs] == [ResponseMatrix(1,si).rs]  ;
            
            perf_all_per_em_per_vs(si,ci,ni) = 100*sum(bhv & ind_h)/sum(ind_h);
            rt_all_per_em_per_vs(si,ci,ni)   = nanmean(ResponseMatrix(1,si).rt(ind_h));
        
            
            
        end
    end
end




%% Plot psychometric function for each emotion on stimulus visual signal

visual_signal = abs(unique([ResponseMatrix(1,:).noise]) - 100);


figure(3);
legends            = {'joy','sadness'};
for pm = 1:ci
    newcolors =[0    0.4470    0.7410
    0.8500    0.3250    0.0980];


    xlabel('Visual Signal (%)')
    ylabel('Mean Performance on all subjects (%)')

     plot(visual_signal, nanmean(squeeze(perf_all_per_em_per_vs(:,pm,:))),'LineWidth',4)
          hold on
      e=errorbar(visual_signal,nanmean(squeeze(perf_all_per_em_per_vs(:,pm,:))),nanstd(squeeze(perf_all_per_em_per_vs(:,pm,:)))/sqrt(size(perf_all_per_em_per_vs,1)),'HandleVisibility','off','Color',newcolors(pm,:),'LineWidth',1);
 
     legend(legends,'Location','northeast')
    xlim([0 100]);
    ylim([0 100]);
    title('Performance')
        set (gca, 'FontSize' ,15)
end
saveas( figure(3) , 'fig3.jpg');
figure(4);
for pm=1:ci

    xlabel('Visual Signal (%)')
    ylabel('Mean Reaction Time on all subjects (sec)')
     plot(visual_signal, nanmean(squeeze(rt_all_per_em_per_vs(:,pm,:))),'LineWidth',4)
        hold on

     errorbar(visual_signal,nanmean(squeeze(rt_all_per_em_per_vs(:,pm,:))),nanstd(squeeze(rt_all_per_em_per_vs(:,pm,:)))/sqrt(size(rt_all_per_em_per_vs,1)),'HandleVisibility','off','Color',newcolors(pm,:),'LineWidth',1)

    legend(legends,'Location','northeast')

    xlim([0 100]);
    title('Reaction Time')
    set (gca, 'FontSize' ,15)
end
saveas( figure(4) , 'fig4.jpg');




%% Sample Fitting logistic regression to sample data

si = 3;
ci = 1;
p_ = round(squeeze(perf_all_per_em_per_tim(si,ci,:))/100*10);

PF = [];
x = times;
y = p_;
n = 10*ones(4,1) ;

b = glmfit(times,[y n],'binomial','link','logit');

yfit = glmval(b,x,'logit','size',n);
figure(1)
plot(x,100* y./n,'ko',x,100*yfit./n,'r-','LineWidth',2)
hold on
plot(x,100*y./n,'LineWidth',1,'Color', 'k')
xlabel('Presentation Time(s)')
ylabel('Percent Correct(%)')
set(gca, 'fontsize', 20);

xlabel('Display Duration (s)');
ylabel('Correct Response(%)');

title('Performance');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
legend('Measured','Fitted')
saveas( figure(1) , 'sample.jpg');
%% Fitting logistic regression to data ( time )
model_p =[];
for si = 1 : size(ResponseMatrix,2)
    for ci = 1: length(emotions)
        p_ = round(squeeze(perf_all_per_em_per_tim(si,ci,:))/100*10);
        
        PF = [];
        x = times;
        y = p_;
        n = 10*ones(4,1) ;
        
        b = glmfit(times,[y n],'binomial','link','probit');
        model_p(si,ci,:) =b;
    
    end
    
end
tresh =-1*squeeze(model_p(:,:,1))./squeeze(model_p(:,:,2));
sensitivity = squeeze(model_p(:,:,2));

%% Removoing out of range data (time)
tresh (tresh >100) =nan; tresh (tresh <-100) =nan;

%% Preview threshold and sestivity of model fitting displayed emotion (time)
figure(1);
 ax = subplot(1,1,1);
hold on

mBar= bar(1:4,squeeze(tresh(:,:)),'FaceColor','g');
ax.XTick = [1:4];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;
xlabel('Emotions')
ylabel('Threshold (s)')
title('Threshold')
set(gca, 'fontsize', 20);


set(gcf,'inverthardcopy','off'); 
saveas( figure(1) , 'threshold-1.jpg');
figure(2);
 ax = subplot(1,1,1);
hold on
bar(1:4,squeeze(sensitivity(:,:)),'FaceColor','g')
ax.XTick = [1:4];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;

xlabel('Emotions')
ylabel('Sensitivity (a.u.)')
title('Sensitivity')
set(gca, 'fontsize', 20);
saveas( figure(2) , 'sensivity-2.jpg');

set(gcf, 'Color', [0.97 0.97 0.97]);



%% Fitting logistic regression to data ( noise )
model_p =[];
for si = 1 : size(ResponseMatrix,2)
    for ci = 1: length(emotions)
        p_ = round(squeeze(perf_all_per_em_per_vs(si,ci,:))/100*10);
        
        PF = [];
        x = visual_signal;
        y = p_;
        n = 10*ones(4,1) ;
        
        b = glmfit(visual_signal,[y n],'binomial','link','probit');
        model_p(si,ci,:) =b;
     
    end
    
end
tresh =-1*squeeze(model_p(:,:,1))./squeeze(model_p(:,:,2));
sensitivity = squeeze(model_p(:,:,2));

%% Removoing out of range data (noise)
tresh (tresh >100) =nan; tresh (tresh <-100) =nan;

%% Preview threshold and sestivity of model fitting displayed emotion (noise)
figure(3);
 ax = subplot(1,1,1);
hold on

mBar= bar(1:4,squeeze(tresh(:,:)),'FaceColor','g');
ax.XTick = [1:4];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;
xlabel('Emotions')
ylabel('Threshold (s)')
title('Threshold')
set(gca, 'fontsize', 20);


set(gcf,'inverthardcopy','off'); 
saveas( figure(3) , 'threshold-3.jpg');




figure(4);
 ax = subplot(1,1,1);
hold on
bar(1:4,squeeze(sensitivity(:,:)),'FaceColor','g')
ax.XTick = [1:4];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;

xlabel('Emotions')
ylabel('Sensitivity (a.u.)')
title('Sensitivity')
set(gca, 'fontsize', 20);
saveas( figure(4) , 'sensivity-4.jpg');

set(gcf, 'Color', [0.97 0.97 0.97]);


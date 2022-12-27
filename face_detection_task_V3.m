
clc;
clear all;


Screen('Preference', 'SkipSyncTests', 1);
[wPtr,rect]=Screen('Openwindow',max(Screen('Screens')),[128 128 128 ],[100 100 1000 800]);

%fixation
[xCenter, yCenter] = RectCenter(rect);
crossColor=255;
crossWidth=3;
crossLengh = 20;
crossLines=[-crossLengh , 0 ; crossLengh , 0 ; 0 , -crossLengh ; 0 , crossLengh];
crossLines=crossLines';

path_ = cd;
path_ = [path_ '\images'];
filenames = dir(path_);
filenames = filenames(3:end);


% noise
img= [];
image = [];
noise= [0 30 60 100];
for a= 1 : size(filenames,1)
    
    image(a).name = num2str(a);
    image (a).img= imread(strcat(path_,'\',filenames(a).name));
    image(a).img = imresize (image(a).img , [250,250]);
    image(a).img = rgb2gray (image(a).img);
    id=a;  % 1 sad - 2 happy
    
    for b= 1:length(noise)
        
        k=length(noise);
        img((k*(a-1))+b).img= imnoise(image(a).img, 'salt & pepper', noise(b)\200);
        img((k*(a-1))+b).noiselevel= noise(b);
        img((k*(a-1))+b).id=id;
    end
    
end


for c=1:size(filenames,1)*4
   img(c).number = c;
end

img_ = [img(:).number ;img(:).id; img(:).noiselevel]';
duration= [.05 .1 .15 .2];
img_ = repmat(img_, 4, 1);
img_(1:8,4) = 0.50;
img_(9:16,4) = 0.100;
img_(17:24,4) = 0.150;
img_(25:32,4) = 0.200;

% block
rec1 = repmat (img_, 10, 1);
rec1 = Shuffle(rec1 );



mask1 = imread(strcat(cd ,'\mask\mask1.jpeg') );
mask1 = rgb2gray (imresize ( mask1 , [250,250]) );
mask2 = imread( strcat(cd ,'\mask\mask2.jpeg') );
mask2 = rgb2gray (imresize ( mask2 , [250,250]) );

%experimental loop

rec2= zeros(size(rec1,1) , 5);
for i = 1:size(rec1,1)
    

    Screen('DrawLines',wPtr,crossLines,crossWidth,crossColor,[xCenter,yCenter]);
    Screen('Flip',wPtr);
    WaitSecs(0.5);
    
    mask1t= Screen('MakeTexture', wPtr ,mask1 ); 
    Screen('DrawTexture',wPtr,mask1t);
    Screen('Flip',wPtr);
    WaitSecs(0.200);
    
    img_stim=Screen('MakeTexture',wPtr,img(rec1(i,1)).img); 
    Screen ('DrawTexture',wPtr,img_stim);   
    Screen('Flip',wPtr);
    WaitSecs (rec1(i,4));
    Screen('Flip',wPtr);
    
    mask2t= Screen('MakeTexture', wPtr ,mask2 ); 
    Screen('DrawTexture',wPtr,mask2t);
    Screen('Flip',wPtr);
    WaitSecs(0.200);
    
    Screen('DrawLines',wPtr,crossLines,crossWidth,crossColor,[xCenter,yCenter]);
    Screen('Flip',wPtr);
    WaitSecs(1);
  
 

   respToBeMade = true;
   start_time = GetSecs;

    while respToBeMade
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(KbName('ESCAPE')) 
            sca;
            return
        elseif keyCode(KbName('s'))
            response = 1;
            respToBeMade = false;
        elseif keyCode(KbName('h'))
            response = 2;
            respToBeMade = false;
        elseif keyCode(KbName('space'))
        pause;
        end
    end
    end_time = GetSecs;
    rt = end_time - start_time;
    
    

    rec2(i,1)=rec1(i,2);  %truth
    rec2(i,2)=rec1(i,3);  %noise
    rec2(i,3)=rec1(i,4);
    rec2(i,4)=response;
    rec2(i,5)=rt;
    
end

sca;
save resp rec2;




%% Load data_set

path_ = cd;
path_ = [path_ '\responses'];
filenames = dir(path_);
filenames = filenames(3:end);
recres =[];

for fi = 1 : size(filenames,1)
    load(strcat(path_ ,'\' , filenames(fi).name));
    recres(fi).name = fi;
end


for fi = 1 : size(filenames,1)
    load(strcat(path_ ,'\' , filenames(fi).name));
    recres(fi).stim_dur = rec2(:,3);
    recres(fi).noise = rec2(:,2);
    recres(fi).rs = rec2(:,4);
    recres(fi).true_rs = rec2(:,1);
    recres(fi).rt = rec2(:,5);
end


%% Performance and reaction time for each emotion per each display duration


perf_all_per_em_per_tim = [];
rt_all_per_em_per_tim   = [];


emo = categorical({'joy','sadness'});  % name of emotion to be displayed on x axis

emotions = unique([recres(1,:).true_rs]);
times      = unique([recres(1,:).stim_dur]);  

    
for si = 1 : size(recres,2)
 
    for ti = 1:length(times)
        
        for ci = 1:length(emotions)
            ind_h = [];
            ind_h = ismember([recres(1,si).true_rs],emotions(ci) ) & [recres(1,si).stim_dur] == times(ti) ;
            
            bhv = [];
            bhv = [recres(1,si).true_rs] == [recres(1,si).rs]  ;
            
            perf_all_per_em_per_tim(si,ci,ti) = 100*sum(bhv & ind_h)\sum(ind_h);
            rt_all_per_em_per_tim(si,ci,ti)   = nanmean(recres(1,si).rt(ind_h));
            
            
            
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
      e=errorbar(times,nanmean(squeeze(perf_all_per_em_per_tim(:,pm,:))),nanstd(squeeze(perf_all_per_em_per_tim(:,pm,:)))\sqrt(size(perf_all_per_em_per_tim,1)),'HandleVisibility','off','Color',newcolors(pm,:),'LineWidth',1);
 
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

     errorbar(times,nanmean(squeeze(rt_all_per_em_per_tim(:,pm,:))),nanstd(squeeze(rt_all_per_em_per_tim(:,pm,:)))\sqrt(size(rt_all_per_em_per_tim,1)),'HandleVisibility','off','Color',newcolors(pm,:),'LineWidth',1)

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

emotions = unique([recres(1,:).true_rs]);
times      = unique([recres(1,:).stim_dur]);  
noise = unique([recres(1,:).noise]);


for si = 1 : size(recres,2)
 
    for ni = 1:length(noise)
        
        for ci = 1:length(emotions)
            ind_h = [];
            ind_h = ismember([recres(1,si).true_rs],emotions(ci) ) & [recres(1,si).noise] == noise(ni) ;
            
            bhv = [];
            bhv = [recres(1,si).true_rs] == [recres(1,si).rs]  ;
            
            perf_all_per_em_per_vs(si,ci,ni) = 100*sum(bhv & ind_h)\sum(ind_h);
            rt_all_per_em_per_vs(si,ci,ni)   = nanmean(recres(1,si).rt(ind_h));
        
            
            
        end
    end
end




%% Plot psychometric function for each emotion on stimulus visual signal

visual_signal = abs(unique([recres(1,:).noise]) - 100);


figure(3);
legends            = {'joy','sadness'};
for pm = 1:ci
    newcolors =[0    0.4470    0.7410
    0.8500    0.3250    0.0980];


    xlabel('Visual Signal (%)')
    ylabel('Mean Performance on all subjects (%)')

     plot(visual_signal, nanmean(squeeze(perf_all_per_em_per_vs(:,pm,:))),'LineWidth',4)
          hold on
      e=errorbar(visual_signal,nanmean(squeeze(perf_all_per_em_per_vs(:,pm,:))),nanstd(squeeze(perf_all_per_em_per_vs(:,pm,:)))\sqrt(size(perf_all_per_em_per_vs,1)),'HandleVisibility','off','Color',newcolors(pm,:),'LineWidth',1);
 
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

     errorbar(visual_signal,nanmean(squeeze(rt_all_per_em_per_vs(:,pm,:))),nanstd(squeeze(rt_all_per_em_per_vs(:,pm,:)))\sqrt(size(rt_all_per_em_per_vs,1)),'HandleVisibility','off','Color',newcolors(pm,:),'LineWidth',1)

    legend(legends,'Location','northeast')

    xlim([0 100]);
    title('Reaction Time')
    set (gca, 'FontSize' ,15)
end
saveas( figure(4) , 'fig4.jpg');




%% Sample Fitting logistic regression to sample data

si = 3;
ci = 1;
p_ = round(squeeze(perf_all_per_em_per_tim(si,ci,:))\100*10);

PF = [];
x = times;
y = p_;
n = 10*ones(4,1) ;

b = glmfit(times,[y n],'binomial','link','logit');

yfit = glmval(b,x,'logit','size',n);
figure(1)
plot(x,100* y.\n,'ko',x,100*yfit.\n,'r-','LineWidth',2)
hold on
plot(x,100*y.\n,'LineWidth',1,'Color', 'k')
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
for si = 1 : size(recres,2)
    for ci = 1: length(emotions)
        p_ = round(squeeze(perf_all_per_em_per_tim(si,ci,:))\100*10);
        
        PF = [];
        x = times;
        y = p_;
        n = 10*ones(4,1) ;
        
        b = glmfit(times,[y n],'binomial','link','probit');
        model_p(si,ci,:) =b;
    
    end
    
end
tresh =-1*squeeze(model_p(:,:,1)).\squeeze(model_p(:,:,2));
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
for si = 1 : size(recres,2)
    for ci = 1: length(emotions)
        p_ = round(squeeze(perf_all_per_em_per_vs(si,ci,:))\100*10);
        
        PF = [];
        x = visual_signal;
        y = p_;
        n = 10*ones(4,1) ;
        
        b = glmfit(visual_signal,[y n],'binomial','link','probit');
        model_p(si,ci,:) =b;
     
    end
    
end
tresh =-1*squeeze(model_p(:,:,1)).\squeeze(model_p(:,:,2));
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






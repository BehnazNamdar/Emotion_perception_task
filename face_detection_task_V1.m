%% This script was written for final exam of Cognition Lab course 


% Task is recognizing emotions with different stimulus duration times and
% also with different noise rates

    % Summary of tasl parameters:
    % ---------------------------

            % fixation time = 500 ms 
            % scrambled image = 200 ms
            % image ( 2 image)
                %interval = 50 , 100 , 150 , 200
                % noise = 0 , 30 , 60 , 100
            % scrambled image = 200 ms
            % delay time = 1 s 
            % response : s & h keys
            % block num = 10 
            % subject num = 4
            
    % What to do with data-sets:
    % --------------------------

            % psychmetric function for each image kind per each visual signal (happy, sad) 
            % psychmetric function for each image kind per each stimulus time (happy, sad)  
            % fit logit function to all the psychmetric function 

    % send "script, mat file, videoscreen , plot " 

%% Sections Outline  


%     1- This script was written for final exam of Cognition Lab course
%     2- Sections Outline 
%     3- preliminary settings
    
    
%     4- make noisy images, save them 
%     5- initiate experiment array
%     6- read scrambled images
%     7- instruction
%     8- experimental loop
     
     
     % READ DATA AND PLOT PSYCHOMETRIC FUNCTIONS
     
%     9- Load data_set
%     10- Performance and reaction time for each emotion per each duration duration
%     11- Plot psychometric function for each emotion on stimulus duration time
%     12- Performance and reaction time for each emotion per each visual signal
%     13- Plot psychometric function for each emotion on stimulus visual signal

     % FITT LOGISTIC REGRESSION TO DATA 
     
%     14- Sample Fitting logistic regression to sample data     
%     15- 
%           15-1- Fitting logistic regression to data - performance ( time )
%           15-2- Removoing out of range data - performance (time)
%           15-3- Preview threshold and sestivity of model fitting displayed emotion - performance (time)
     
%     16-
%           16-1- Fitting logistic regression to data - rt ( time )
%           16-2- Removoing out of range data - rt (time)
%           16-3- Preview threshold and sestivity of model fitting displayed emotion - rt (time)
     
%     17-
%           17-1- Fitting logistic regression to data - performance ( noise )
%           17-2- Removoing out of range data - performance (noise)
%           17-3- Preview threshold and sestivity of model fitting displayed emotion - performance (noise)
     
%     18-
%           18-1- Fitting logistic regression to data - rt ( noise )
%           18-2 Removoing out of range data - rt (noise)
%           18-3- Preview threshold and sestivity of model fitting displayed emotion - rt (noise)
    
%%  preliminary settings 
    
% Clear the workspace
clear all;close all;sca;clc;


% Enter subject code or number
sub_code = input ('Enter subject code or number');


% -------------------------
% Initialize screens
% -------------------------

Screen('Preference', 'SkipSyncTests', 1);
[width, height]=Screen('WindowSize',max(Screen('Screens')));
rect = [0 0 width height];

screenNumber = max(Screen('Screens'));
[wPtr, windowRect] = Screen('OpenWindow', screenNumber, 128, rect );



% white=WhiteIndex(screenNumber);
% black=BlackIndex(screenNumber);


% -------------------------
% set fixation point parameters
% -------------------------

[xCenter, yCenter] = RectCenter(windowRect);
crossColor=255;
crossWidth=3;
crossLengh = 20;
crossLines=[-crossLengh , 0 ; crossLengh , 0 ; 0 , -crossLengh ; 0 , crossLengh];
crossLines=crossLines';




% ---------------------
% Keyboard information
% ---------------------

%Keypress settings
KbName('UnifyKeyNames');
KbDeviceIndex = [];

PauseKey = KbName('space');
escapeKey = KbName('ESCAPE');
sad_Key = KbName('s');
happy_Key = KbName('h' );

HideCursor;


% -------------------------
% set stimulus parameters
% -------------------------

% image path
path_ = cd;
path_ = [path_ '/images'];
files = dir (path_);
files = files(3: end);


fixation_time_dur = 0.5  ; % s
scramb_image_time_dur = 0.2; % s
img_noise = [0 , 30 , 60 , 100 ]; %
img_noise_num = length(img_noise); 
stimulus_time_dur = [0.05 , 0.1 , 0.15 , 0.2]; % s
stimulus_time_dur_num = length(stimulus_time_dur);
delay_time_dur = 1 ; % s 
block_num = 10 ;
trial_num_in_block = size(files,1) .* img_noise_num .* stimulus_time_dur_num  ;
total_trial = trial_num_in_block * block_num;
stimulus_size = [400 400];




% make new dir for new images
mkdir new_images;
newdir = [cd , '/new_images'] ;
newfiles = dir (newdir);
newfiles = newfiles(3: end);



%% make noisy images, save them 
k = 1; 
for i = 1 : size(files,1)
    image = imread (strcat ( path_ , '/' , files(i).name));  % image read
    image = rgb2gray (imresize ( image , stimulus_size) ); % make image graysclale and resize it
    
    
    for j = 1 : img_noise_num
        
        image = imnoise (image ,'salt & pepper', 0.01 * img_noise(j) );
        image_name = [newdir,'/Image_',num2str(i),'_', num2str(j) ,'.jpg'] ;
        imwrite(image , image_name );
        
        res(k).noise_level = img_noise(j); 
        if i == 1
            res(k).imtype = sad_Key;
        elseif i == 2 
            res(k).imtype = happy_Key;
        end
        res(k).pic = image;
        k = k + 1;
    end
    
end
i = 0;
j = 0;
k = 0;

%% initiate experiment array 


a = repmat( [res] ,stimulus_time_dur_num , 1);

for i = 1:stimulus_time_dur_num
    for j = 1: length(a)
        a(i , j).stimulus_time_dur = stimulus_time_dur(i);
    end
end
i = 0; j = 0;
a = reshape(a , [trial_num_in_block,1]);

%% read scrambled images 

scramb_1 = imread( strcat(cd , '/' , 'imgMask/' , '‌B_01.jpeg') );
scramb_1 = rgb2gray (imresize ( scramb_1 , stimulus_size) );

scramb_2 = imread( strcat(cd , '/' , 'imgMask/' , 'B_02.jpeg') );
scramb_2 = rgb2gray (imresize ( scramb_2 , stimulus_size) );


%% instruction

Screen('TextFont', wPtr, 'Times'); % set the font for the screen to Times
Screen('TextSize', wPtr, 30); % set the font size for the screen to 30
text = 'You will see numbers of images containing either sad or happy faces.\n You should press "s" key when you see sad faces and \n press "h" key when you see happy faces.\n Press any key to start when you are ready';
 

DrawFormattedText(wPtr, text ,'center', 'center', 255);
Screen('Flip', wPtr);
KbWait;

%% experimental loop

blocks = 1;
while blocks <= block_num
    %shuffle the basic structure to start a trial
    a = Shuffle(a);
    for k = 1 : trial_num_in_block 

        % Draw fixation point

        Screen('DrawLines',wPtr,crossLines,crossWidth,crossColor,[xCenter,yCenter]);
        Screen('Flip',wPtr);
        WaitSecs(fixation_time_dur);


        % show scrambled image1 

        scramb_1_ = Screen('MakeTexture', wPtr , scramb_1 ); 
        Screen('DrawTexture',wPtr,scramb_1_);
        Screen('Flip',wPtr);
        WaitSecs(scramb_image_time_dur);



        % show stimulus 

        stim = Screen('MakeTexture', wPtr , a(k).pic); 
        Screen('DrawTexture',wPtr,stim);
        Screen('Flip',wPtr);
        WaitSecs(a(k).stimulus_time_dur);


        % show scrambled image2 

        scramb_2_ = Screen('MakeTexture', wPtr , scramb_2 ); 
        Screen('DrawTexture',wPtr,scramb_2_);
        Screen('Flip',wPtr);
        WaitSecs(scramb_image_time_dur);

        % delay 

        Screen('DrawLines',wPtr,crossLines,crossWidth,crossColor,[xCenter,yCenter]);
        Screen('Flip',wPtr);
        WaitSecs(delay_time_dur);
        Screen('Flip',wPtr);



        RestrictKeysForKbCheck([ escapeKey, happy_Key, sad_Key , PauseKey] );

        keyIsDown = 0;
        tStart = GetSecs;
        while ~keyIsDown 
            [keyIsDown,secs, keyCode] = KbCheck(-1);
        end

        pressedKey = find(keyCode);
            if keyCode(escapeKey)
                ShowCursor;
                sca;
                return;
            elseif keyCode(PauseKey)
                pause;
            elseif keyCode(happy_Key)
                response = 11;
                keyIsDown = 0;
            elseif keyCode(sad_Key)
                response = 22;
                keyIsDown = 0;
            end

        % calculate reaction time 
        rt = secs - tStart;
       
        % response array
        tdata(((blocks-1)*trial_num_in_block) + k).rt = rt;
        tdata(((blocks-1)*trial_num_in_block) + k).rs = response;
        tdata(((blocks-1)*trial_num_in_block) + k).noise_level = a(k).noise_level;
        tdata(((blocks-1)*trial_num_in_block) + k).true_rs = a(k).imtype;
        tdata(((blocks-1)*trial_num_in_block) + k).stim_dur = a(k).stimulus_time_dur;
        tdata(((blocks-1)*trial_num_in_block) + k).block = blocks;

    end
    blocks = blocks + 1;
end
    

save mat_file tdata;
ShowCursor;
sca;
close all; clc; clearvars;


%% Load data_set

path_ = cd;
path_ = [path_ '/Final_exam_01'];
filenames = dir(path_);
filenames = filenames(3:end);
data_set =[];

for fi = 1 : size(filenames,1)
    load(strcat(path_ ,'/' , filenames(fi).name));
    data_set(fi).subject = extract(filenames(fi).name, digitsPattern );
    data_set(fi).stim_dur = ([tdata.stim_dur]);
    data_set(fi).noise = [tdata.noise_level];
    data_set(fi).rs = [tdata.rs];
    data_set(fi).true_rs = [tdata.true_rs];
    data_set(fi).rt = [tdata.rt];
    data_set(fi).block = [tdata.block];
end




%% Performance and reaction time for each emotion per each display duration


perf_all_per_em_per_tim = [];
rt_all_per_em_per_tim   = [];
rt_all_per_em_per_tim_cr =[];
rt_all_per_em_per_tim_wr =[];
cl={};

emo = categorical({'joy','sadness'});  % name of emotion to be displayed on x axis

emotions = unique([data_set(1,:).true_rs]);
times      = unique([data_set(1,:).stim_dur]);  

    
for si = 1 : size(data_set,2)
 
    for ti = 1:length(times)
        
        for ci = 1:length(emotions)
            ind_h = [];
            ind_h = ismember([data_set(1,si).true_rs],emotions(ci) ) & [data_set(1,si).stim_dur] == times(ti) ;
            
            bhv = [];
            bhv = [data_set(1,si).true_rs] == [data_set(1,si).rs]  ;
            
            perf_all_per_em_per_tim(si,ci,ti) = 100*sum(bhv & ind_h)/sum(ind_h);
            rt_all_per_em_per_tim(si,ci,ti)   = nanmean(data_set(1,si).rt(ind_h));
            rt_all_per_em_per_tim_cr(si,ci,ti)   = nanmean(data_set(1,si).rt(ind_h & bhv));
            rt_all_per_em_per_tim_wr(si,ci,ti)   = nanmean(data_set(1,si).rt(ind_h & ~bhv));
            
            
        end
    end
end


%% Plot psychometric function for each emotion on stimulus duration time

figure(1)
legends            = {'joy','sadness'};
for pm = 1:ci
    newcolors =[0    0.4470    0.7410
    0.8500    0.3250    0.0980];
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840

%     subplot(1,2,1);
    xlabel('Display duration (s)')
    ylabel('Mean Performance on all subjects (%)')

     plot(times, nanmean(squeeze(perf_all_per_em_per_tim(:,pm,:))),'LineWidth',4)
          hold on
      e=errorbar(times,nanmean(squeeze(perf_all_per_em_per_tim(:,pm,:))),nanstd(squeeze(perf_all_per_em_per_tim(:,pm,:)))/sqrt(size(perf_all_per_em_per_tim,1)),'HandleVisibility','off','Color',newcolors(pm,:),'LineWidth',1);
 
     legend(legends,'Location','northwest')
    xlim([0 0.25]);
    ylim([20 100]);
    title('Performance')
        set (gca, 'FontSize' ,15)
end
saveas( figure(1) , 'fig1.jpg');

figure(2)
for pm=1:ci
%     subplot(1,2,2);
    xlabel('Display diuration (s)')
    ylabel('Mean Reaction Time on all subjects (sec)')
     plot(times, nanmean(squeeze(rt_all_per_em_per_tim(:,pm,:))),'LineWidth',4)
        hold on

     errorbar(times,nanmean(squeeze(rt_all_per_em_per_tim(:,pm,:))),nanstd(squeeze(rt_all_per_em_per_tim(:,pm,:)))/sqrt(size(rt_all_per_em_per_tim,1)),'HandleVisibility','off','Color',newcolors(pm,:),'LineWidth',1)

    legend(legends,'Location','northwest')
%     ylim([0,100]);
    xlim([0 0.25]);
    title('Reaction Time')
    set (gca, 'FontSize' ,15)
end
saveas( figure(2) , 'fig2.jpg');


%% Performance and reaction time for each emotion per each visual signal


perf_all_per_em_per_vs = [];
rt_all_per_em_per_vs   = [];
rt_all_per_em_per_vs_cr =[];
rt_all_per_em_per_vs_wr =[];
cl={};

emo = categorical({'joy','sadness'});  % name of emotion to be displayed on x axis

emotions = unique([data_set(1,:).true_rs]);
times      = unique([data_set(1,:).stim_dur]);  
noise = unique([data_set(1,:).noise]);


for si = 1 : size(data_set,2)
 
    for ni = 1:length(noise)
        
        for ci = 1:length(emotions)
            ind_h = [];
            ind_h = ismember([data_set(1,si).true_rs],emotions(ci) ) & [data_set(1,si).noise] == noise(ni) ;
            
            bhv = [];
            bhv = [data_set(1,si).true_rs] == [data_set(1,si).rs]  ;
            
            perf_all_per_em_per_vs(si,ci,ni) = 100*sum(bhv & ind_h)/sum(ind_h);
            rt_all_per_em_per_vs(si,ci,ni)   = nanmean(data_set(1,si).rt(ind_h));
            rt_all_per_em_per_vs_cr(si,ci,ni)   = nanmean(data_set(1,si).rt(ind_h & bhv));
            rt_all_per_em_per_vs_wr(si,ci,ni)   = nanmean(data_set(1,si).rt(ind_h & ~bhv));
            
            
        end
    end
end




%% Plot psychometric function for each emotion on stimulus visual signal

visual_signal = abs(unique([data_set(1,:).noise]) - 100);


figure(3);
legends            = {'joy','sadness'};
for pm = 1:ci
    newcolors =[0    0.4470    0.7410
    0.8500    0.3250    0.0980];
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840

    xlabel('Visual Signal (%)')
    ylabel('Mean Performance on all subjects (%)')

     plot(visual_signal, nanmean(squeeze(perf_all_per_em_per_vs(:,pm,:))),'LineWidth',4)
          hold on
      e=errorbar(visual_signal,nanmean(squeeze(perf_all_per_em_per_vs(:,pm,:))),nanstd(squeeze(perf_all_per_em_per_vs(:,pm,:)))/sqrt(size(perf_all_per_em_per_vs,1)),'HandleVisibility','off','Color',newcolors(pm,:),'LineWidth',1);
 
     legend(legends,'Location','northwest')
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

    legend(legends,'Location','northwest')
%     ylim([0,100]);
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
% ylim([0,1]);
% legend(legends,'Location','southeast')
title('Performance');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
legend('Measured','Fitted')
saveas( figure(1) , 'samplefitting.jpg');

%% Fitting logistic regression to data - performance ( time )
model_p =[];
for si = 1 : size(data_set,2)
    for ci = 1: length(emotions)
        p_ = round(squeeze(perf_all_per_em_per_tim(si,ci,:))/100*10);
        
        PF = [];
        x = times;
        y = p_;
        n = 10*ones(4,1) ;
        
        b = glmfit(times,[y n],'binomial','link','probit');
        model_p(si,ci,:) =b;
        % yfit = glmval(b,x,'probit','size',n);
        % plot(x, y./n,'o',x,yfit./n,'-','LineWidth',2)
    end
    
end
tresh =-1*squeeze(model_p(:,:,1))./squeeze(model_p(:,:,2));
sensitivity = squeeze(model_p(:,:,2));

%% Removoing out of range data - performance (time)
tresh (tresh >100) =nan; tresh (tresh <-100) =nan;

%% Preview threshold and sestivity of model fitting displayed emotion - performance (time)
figure(1);
 ax = subplot(1,1,1);
hold on

mBar= bar(1:2,squeeze(tresh(:,:)),'FaceColor','g');
% errorbar(1:2,nanmean(tresh(:,:)),nanstd(tresh(:,:))/sqrt(size(tresh,1)),'LineWidth',2,'Color','k','LineStyle','none')
ax.XTick = [1:2];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;
xlabel('Emotions')
ylabel('Threshold (s)')
title('Threshold')
set(gca, 'fontsize', 20);
% axis ([0 7 -5 70])
% mBar.BarWidth=0.9;
% mBar.LineWidth=.5
% mBar.FaceColor=[0 0.6 0.6];
% set(gca,'YGrid','on');
% set(gca,'XGrid','on');

set(gcf,'inverthardcopy','off'); 
saveas( figure(1) , 'perf_threshold_time.jpg');
figure(2);
 ax = subplot(1,1,1);
hold on
bar(1:2,squeeze(sensitivity(:,:)),'FaceColor','g')
% errorbar(1:2,nanmean(sensitivity(:,:)),nanstd(sensitivity(:,:))/sqrt(size(sensitivity,1)),'LineWidth',2,'Color','k','LineStyle','none')
ax.XTick = [1:2];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;

xlabel('Emotions')
ylabel('Sensitivity (a.u.)')
title('Sensitivity')
set(gca, 'fontsize', 20);
% axis ([0 7 -0.05 0.5])
saveas( figure(2) , 'perf_sensivity_time.jpg');

set(gcf, 'Color', [0.97 0.97 0.97]);




%% Fitting logistic regression to data - rt ( time )
model_p =[];
for si = 1 : size(data_set,2)
    for ci = 1: length(emotions)
        p_ = round(squeeze(rt_all_per_em_per_tim(si,ci,:))/100*10);
        
        PF = [];
        x = times;
        y = p_;
        n = 10*ones(4,1) ;
        
        b = glmfit(times,[y n],'binomial','link','probit');
        model_p(si,ci,:) =b;
        % yfit = glmval(b,x,'probit','size',n);
        % plot(x, y./n,'o',x,yfit./n,'-','LineWidth',2)
    end
    
end
tresh =-1*squeeze(model_p(:,:,1))./squeeze(model_p(:,:,2));
sensitivity = squeeze(model_p(:,:,2));

%% Removoing out of range data - rt(time)
tresh (tresh >100) =nan; tresh (tresh <-100) =nan;

%% Preview threshold and sestivity of model fitting displayed emotion - rt (time)
figure(1);
 ax = subplot(1,1,1);
hold on

mBar= bar(1:2,squeeze(tresh(:,:)),'FaceColor','g');
% errorbar(1:2,nanmean(tresh(:,:)),nanstd(tresh(:,:))/sqrt(size(tresh,1)),'LineWidth',2,'Color','k','LineStyle','none')
ax.XTick = [1:2];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;
xlabel('Emotions')
ylabel('Threshold (s)')
title('Threshold')
set(gca, 'fontsize', 20);
% axis ([0 7 -5 70])
% mBar.BarWidth=0.9;
% mBar.LineWidth=.5
% mBar.FaceColor=[0 0.6 0.6];
% set(gca,'YGrid','on');
% set(gca,'XGrid','on');

set(gcf,'inverthardcopy','off'); 
saveas( figure(1) , 'rt_threshold_time.jpg');
figure(2);
 ax = subplot(1,1,1);
hold on
bar(1:2,squeeze(sensitivity(:,:)),'FaceColor','g')
% errorbar(1:2,nanmean(sensitivity(:,:)),nanstd(sensitivity(:,:))/sqrt(size(sensitivity,1)),'LineWidth',2,'Color','k','LineStyle','none')
ax.XTick = [1:2];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;

xlabel('Emotions')
ylabel('Sensitivity (a.u.)')
title('Sensitivity')
set(gca, 'fontsize', 20);
% axis ([0 7 -0.05 0.5])
saveas( figure(2) , 'rt_sensivity_time.jpg');

set(gcf, 'Color', [0.97 0.97 0.97]);



%% Fitting logistic regression to data - performance ( noise )
model_p =[];
for si = 1 : size(data_set,2)
    for ci = 1: length(emotions)
        p_ = round(squeeze(perf_all_per_em_per_vs(si,ci,:))/100*10);
        
        PF = [];
        x = visual_signal;
        y = p_;
        n = 10*ones(4,1) ;
        
        b = glmfit(visual_signal,[y n],'binomial','link','probit');
        model_p(si,ci,:) =b;
        % yfit = glmval(b,x,'probit','size',n);
        % plot(x, y./n,'o',x,yfit./n,'-','LineWidth',2)
    end
    
end
tresh =-1*squeeze(model_p(:,:,1))./squeeze(model_p(:,:,2));
sensitivity = squeeze(model_p(:,:,2));

%% Removoing out of range data - performance (noise)
tresh (tresh >100) =nan; tresh (tresh <-100) =nan;

%% Preview threshold and sestivity of model fitting displayed emotion - performance (noise)
figure(3);
 ax = subplot(1,1,1);
hold on

mBar= bar(1:2,squeeze(tresh(:,:)),'FaceColor','g');
% errorbar(1:2,nanmean(tresh(:,:)),nanstd(tresh(:,:))/sqrt(size(tresh,1)),'LineWidth',2,'Color','k','LineStyle','none')
ax.XTick = [1:2];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;
xlabel('Emotions')
ylabel('Threshold (s)')
title('Threshold')
set(gca, 'fontsize', 20);
% axis ([0 7 -5 70])
% mBar.BarWidth=0.9;
% mBar.LineWidth=.5
% mBar.FaceColor=[0 0.6 0.6];
% set(gca,'YGrid','on');
% set(gca,'XGrid','on');

set(gcf,'inverthardcopy','off'); 
saveas( figure(3) , 'perf_threshold_noise.jpg');
figure(4);
 ax = subplot(1,1,1);
hold on
bar(1:2,squeeze(sensitivity(:,:)),'FaceColor','g')
% errorbar(1:2,nanmean(sensitivity(:,:)),nanstd(sensitivity(:,:))/sqrt(size(sensitivity,1)),'LineWidth',2,'Color','k','LineStyle','none')
ax.XTick = [1:2];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;

xlabel('Emotions')
ylabel('Sensitivity (a.u.)')
title('Sensitivity')
set(gca, 'fontsize', 20);
% axis ([0 7 -0.05 0.5])
saveas( figure(4) , 'perf_sensivity_noise.jpg');

set(gcf, 'Color', [0.97 0.97 0.97]);



%% Fitting logistic regression to data - rt( noise )
model_p =[];
for si = 1 : size(data_set,2)
    for ci = 1: length(emotions)
        p_ = round(squeeze(rt_all_per_em_per_vs(si,ci,:))/100*10);
        
        PF = [];
        x = visual_signal;
        y = p_;
        n = 10*ones(4,1) ;
        
        b = glmfit(visual_signal,[y n],'binomial','link','probit');
        model_p(si,ci,:) =b;
        % yfit = glmval(b,x,'probit','size',n);
        % plot(x, y./n,'o',x,yfit./n,'-','LineWidth',2)
    end
    
end
tresh =-1*squeeze(model_p(:,:,1))./squeeze(model_p(:,:,2));
sensitivity = squeeze(model_p(:,:,2));

%% Removoing out of range data - rt (noise)
tresh (tresh >100) =nan; tresh (tresh <-100) =nan;

%% Preview threshold and sestivity of model fitting displayed emotion -rt (noise)
figure(3);
 ax = subplot(1,1,1);
hold on

mBar= bar(1:2,squeeze(tresh(:,:)),'FaceColor','g');
% errorbar(1:2,nanmean(tresh(:,:)),nanstd(tresh(:,:))/sqrt(size(tresh,1)),'LineWidth',2,'Color','k','LineStyle','none')
ax.XTick = [1:2];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;
xlabel('Emotions')
ylabel('Threshold (s)')
title('Threshold')
set(gca, 'fontsize', 20);
% axis ([0 7 -5 70])
% mBar.BarWidth=0.9;
% mBar.LineWidth=.5
% mBar.FaceColor=[0 0.6 0.6];
% set(gca,'YGrid','on');
% set(gca,'XGrid','on');

set(gcf,'inverthardcopy','off'); 
saveas( figure(3) , 'rt_threshold_noise.jpg');
figure(4);
 ax = subplot(1,1,1);
hold on
bar(1:2,squeeze(sensitivity(:,:)),'FaceColor','g')
% errorbar(1:2,nanmean(sensitivity(:,:)),nanstd(sensitivity(:,:))/sqrt(size(sensitivity,1)),'LineWidth',2,'Color','k','LineStyle','none')
ax.XTick = [1:2];
ax.XTickLabel = legends;
ax.XTickLabelRotation = 45;

xlabel('Emotions')
ylabel('Sensitivity (a.u.)')
title('Sensitivity')
set(gca, 'fontsize', 20);
% axis ([0 7 -0.05 0.5])
saveas( figure(4) , 'rt_sensivity_noise.jpg');

set(gcf, 'Color', [0.97 0.97 0.97]);
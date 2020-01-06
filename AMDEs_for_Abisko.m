%% Villum Research Station
% load Hg data 
load('Hg_VRS_all.mat')

% Defining AMDE as Hg less than 0.5 and three consecutive decreasing measurements
AMDE_VRS=[];
for x=1:length(HG)-3
    if HG(x,2)<0.5 && (HG(x,2)>HG(x+1,2) && HG(x+1,2)>HG(x+2,2))
        AMDE_VRS=[AMDE_VRS;floor(HG(x,1))];
    end
end
% Individual Days and assigning years
AMDE_VRS=unique(AMDE_VRS);
AMDE_VRS(:,2)=year(datetime(AMDE_VRS(:,1),'ConvertFrom','datenum'));
AMDE_VRS(:,3)=month(datetime(AMDE_VRS(:,1),'ConvertFrom','datenum'));

% Filtering for only 2011-2014
AMDE_VRS_11_14=[];
for x=2011:2014
    AMDE_VRS_11_14=[AMDE_VRS_11_14;AMDE_VRS(AMDE_VRS(:,2)==x,:)];
end

%% Ny Ålesund Data
% Loading and Cleaning Data
Hg_NYA_11=load('Hg_NYA_2011_12.txt');
Hg_NYA_11(:,1)=doy2date(Hg_NYA_11(:,1),ones(length(Hg_NYA_11),1)*2011)+1;

Hg_NYA_12=load('Hg_NYA_2012_13.txt');
Hg_NYA_12(:,1)=doy2date(Hg_NYA_12(:,1),ones(length(Hg_NYA_12),1)*2012)+1;

Hg_NYA_13=load('Hg_NYA_2013_14.txt');
Hg_NYA_13(:,1)=doy2date(Hg_NYA_13(:,1),ones(length(Hg_NYA_13),1)*2013)+1;

Hg_NYA_14=load('Hg_NYA_2014_15.txt');
Hg_NYA_14(:,1)=doy2date(Hg_NYA_14(:,1),ones(length(Hg_NYA_14),1)*2014)+1;

Hg_NYA_all=[Hg_NYA_11;Hg_NYA_12;Hg_NYA_13;Hg_NYA_14;];
Hg_NYA_all(:,2)=[];
Hg_NYA_all(:,3)=[];
Hg_NYA_all(Hg_NYA_all(:,2)==99.9990,2)=NaN;
%% Atmospheric Mercury Depletion Events
% Defining AMDE as Hg less than 0.5 and twwo consecutive decreasing measurements
AMDE_NYA11_14=[];
for x=1:length(Hg_NYA_all)-3
    if Hg_NYA_all(x,2)<0.5 && (Hg_NYA_all(x,2)>Hg_NYA_all(x+1,2))
        AMDE_NYA11_14=[AMDE_NYA11_14;floor(Hg_NYA_all(x,1))];
    end
end

% Finding unique days and adding year and month for indexing
AMDE_NYA11_14=unique(AMDE_NYA11_14);
AMDE_NYA11_14(:,2)=year(datetime(AMDE_NYA11_14(:,1),'ConvertFrom','datenum'));
AMDE_NYA11_14(:,3)=month(datetime(AMDE_NYA11_14(:,1),'ConvertFrom','datenum'));
%% Stats for VRS
% For calculating occurence frequency
AMDE_VRS_11_14(:,4)=1; 
AMDE_VRS_11_14=array2table(AMDE_VRS_11_14);
AMDE_VRS_11_14.Properties.VariableNames={'Date','Year','Month','AMDE'};
% Group Statistics
Stats_AMDE_VRS=grpstats(AMDE_VRS_11_14,{'Year','Month'},[],'DataVars',{'Year','Month'});

%% Stats for NYA
% For calculating occurence frequency
AMDE_NYA11_14(:,4)=1; 
AMDE_NYA11_14=array2table(AMDE_NYA11_14);
AMDE_NYA11_14.Properties.VariableNames={'Date','Year','Month','AMDE'};
% Group Statistics
Stats_AMDE_NYA=grpstats(AMDE_NYA11_14,{'Year','Month'},[],'DataVars',{'Year','Month'});
%% Putting data into a matrix for plotting with bar plot
AMDE_VRS_Freq=NaN(12,4);
b=[2011;2012;2013;2014];
for a=1:size(Stats_AMDE_VRS,1)
    for x=1:length(b)
        for y=1:12
                if Stats_AMDE_VRS{a,1}==b(x) && Stats_AMDE_VRS{a,2}==y
                    AMDE_VRS_Freq(y,x)=Stats_AMDE_VRS{a,3};
                   
                end
        end
    end
end
%% Putting data into a matrix for plotting with bar plot
AMDE_NYA_Freq=NaN(12,4);
b=[2011;2012;2013;2014];
for a=1:size(Stats_AMDE_NYA,1)
    for x=1:length(b)
        for y=1:12
                if Stats_AMDE_NYA{a,1}==b(x) && Stats_AMDE_NYA{a,2}==y
                    AMDE_NYA_Freq(y,x)=Stats_AMDE_NYA{a,3};
                   
                end
        end
    end
end

%% Plotting
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
bar(AMDE_VRS_Freq,'LineWidth',1.5)
legend({'2011','2012','2013','2014'},'FontSize',26,'NumColumns',2)
xticklabels({'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'});
ylabel('# AMDE Days')
title('Villum Research Station')
ylim([0 10])
ax=gca;
ax.FontSize=18;
ax.FontWeight='bold';
ax.YTick=[0:2:10];

subplot(2,1,2)
bar(AMDE_NYA_Freq,'LineWidth',1.5)
legend({'2011','2012','2013','2014'},'FontSize',26,'NumColumns',2)
xticklabels({'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'});
ylabel('# AMDE Days')
title('Zeppelin Mountain')
ylim([0 10])
ax=gca;
ax.FontSize=18;
ax.FontWeight='bold';
ax.YTick=[0:2:10];
%% Saving data, it is commented out since I already saved it but didn't want to delete it
% AMDE_VRS=table(datetime(AMDE11_14.(1),'ConvertFrom','datenum','Format','yyyy-MM-dd HH:mm:ss'));
% AMDE_NYA=table(datetime(AMDE_NYA11_14.(1),'ConvertFrom','datenum','Format','yyyy-MM-dd HH:mm:ss'));
% 
% writetable(AMDE_VRS,'AMDE_VRS.txt','WriteVariableNames',false)
% writetable(AMDE_NYA,'AMDE_NYA.txt','WriteVariableNames',false)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Plotting each AMDE day for NYA
% Here I plotted the GEM time series for each AMDE, I commented it out
% since it producing a lot of graphs

AMDE_NYA11_14.(1)=datenum(AMDE_NYA11_14.(1));
AMDE_NYA11_14=table2array(AMDE_NYA11_14);
%for x=1:length(AMDE_NYA11_14)
    figure 
    plot(Hg_NYA_all(floor(Hg_NYA_all(:,1))==AMDE_NYA11_14(x),2),'k','LineWidth',1.5);
    xlabel('Hour of Day')
    ylabel('Hg / ng m^{-3}')
    xlim([0 24])
%     xticks([0:4:24])
%     xticklabels([0:4:24])
    title(datestr(AMDE_NYA11_14(x)))
    ax=gca;
    ax.FontWeight='bold';
    ax.FontSize=12;
%end

%%
    figure 
    AMDE_Jan14=Hg_NYA_all(floor(Hg_NYA_all(:,1))==datenum(2014,1,7) | floor(Hg_NYA_all(:,1))==datenum(2014,1,8),:);
    plot(AMDE_Jan14(:,1),AMDE_Jan14(:,2),'k','LineWidth',1.5);
    datetick('x',15)
  
    ylabel('GEM / ng m^{-3}')

   
%     xticks([0:4:24])
%     xticklabels([0:4:24])
    title('Jan 7 & 8, 2014 - Zep')
    ax=gca;
    ax.FontWeight='bold';
    ax.FontSize=12;





%% Plotting each AMDE day for VRS
 
% AMDE_VRS_11_14.(1)=datenum(AMDE_VRS_11_14.(1));
% AMDE_VRS_11_14=table2array(AMDE_VRS_11_14);
% for x=1:length(AMDE_VRS_11_14)
%     figure 
%     plot(HG(floor(HG(:,1))==AMDE_VRS_11_14(x),2),'k','LineWidth',1.5);
%     xlabel('Hour of Day')
%     ylabel('Hg / ng m^{-3}')
%     title(datestr(AMDE_VRS_11_14(x)))
%     ax=gca;
%     ax.FontWeight='bold';
%     ax.FontSize=12;
% end
%% Plotting yearly Hg for VRS
% Formatting data into day of year for daily stats over the years
hg_vrs_11_14=HG;
st=find(hg_vrs_11_14(:,1)==datenum(2011,1,1));
en=find(hg_vrs_11_14(:,1)==datenum(2014,12,31));
hg_vrs_11_14=hg_vrs_11_14(st:en,:);
hg_vrs_11_14_doy=hg_vrs_11_14;
doy=date2doy(hg_vrs_11_14(:,1));
doy=floor(doy);
hg_vrs_11_14_doy(:,1)=doy;
 
% The data was normally distributed so the mean and std are appropriate
% quantities 
yearly_HGmean_VRS=NaN(365,2);
yearly_HGstd_VRS=NaN(365,2);
for x=1:365
    yearly_HGmean_VRS(x,2)=mean(hg_vrs_11_14_doy(hg_vrs_11_14_doy(:,1)==x,2),'omitnan');
    yearly_HGmean_VRS(x,1)=x;
    yearly_HGstd_VRS(x,1)=x;
    yearly_HGstd_VRS(x,2)=std(hg_vrs_11_14_doy(hg_vrs_11_14_doy(:,1)==x,2),'omitnan');
end
jun_aug_hg_vrs=month(datetime(hg_vrs_11_14(:,1),'ConvertFrom','datenum')) >5 & month(datetime(hg_vrs_11_14(:,1),'ConvertFrom','datenum'))<9;
hg_vrs_11_14_jun_aug_mean=mean(hg_vrs_11_14(jun_aug_hg_vrs,2));
hg_vrs_11_14_jun_aug_max=max(hg_vrs_11_14(jun_aug_hg_vrs,2));
% figure
% plot(yearly_HGmean_VRS(:,1),yearly_HGmean_VRS(:,2),'LineWidth',2)
% hold on
% errorbar(yearly_HGmean_VRS(:,1),yearly_HGmean_VRS(:,2),yearly_HGstd_VRS(:,2))

%% Plotting yearly Hg for NYA
hg_NYA=Hg_NYA_all;
doy=date2doy(hg_NYA(:,1));
doy=floor(doy);
hg_NYA_doy=hg_NYA;
hg_NYA_doy(:,1)=doy;
 
yearly_HGmean_NYA=NaN(365,2);
yearly_HGstd_NYA=NaN(365,2);
for x=1:365
    yearly_HGmean_NYA(x,2)=mean(hg_NYA_doy(hg_NYA_doy(:,1)==x,2),'omitnan');
    yearly_HGmean_NYA(x,1)=x;
    yearly_HGstd_NYA(x,1)=x;
    yearly_HGstd_NYA(x,2)=std(hg_NYA_doy(hg_NYA_doy(:,1)==x,2),'omitnan');
end
jun_aug_hg_nya=month(datetime(Hg_NYA_all(:,1),'ConvertFrom','datenum')) >5 & month(datetime(Hg_NYA_all(:,1),'ConvertFrom','datenum'))<9;
hg_nya_11_14_jun_aug_mean=mean(Hg_NYA_all(jun_aug_hg_nya,2),'omitnan');
hg_nya_11_14_jun_aug_max=max(Hg_NYA_all(jun_aug_hg_nya,2));
% figure
% plot(yearly_HGmean_NYA(:,1),yearly_HGmean_NYA(:,2),'k','LineWidth',3)
% hold on
% errorbar(yearly_HGmean_NYA(:,1),yearly_HGmean_NYA(:,2),yearly_HGstd_NYA(:,2));
%% Importing, cleaning, and plotting Sea Ice
addpath('C:\Users\au595124\Documents\PhD Thesis\Courses\Abisko')

SI_NYA=readtable('AMDE_SI_NYA_11_14.csv');
SI_NYA=removevars(SI_NYA,'Var1');

SI_VRS=readtable('AMDE_SI_vrs_11_14.csv');
SI_VRS=removevars(SI_VRS,'Var1');

figure
plot(SI_VRS.time,SI_VRS.SI)
hold on
plot(SI_NYA.time,SI_NYA.SI)
hold off

%% Plotting daily SI median for NYA
% Formatting, averaging, and plotting data for daily mean for the period
SI_NYA.(1)=datenum(SI_NYA.(1));
NYA_SI=table2array(SI_NYA);
doy_nya_si=date2doy(NYA_SI(:,1));
doy_nya_si=floor(doy_nya_si);
NYA_SI(:,1)=doy_nya_si;
 
yearly_SImean_NYA=NaN(365,2);
yearly_SIstd_NYA=NaN(365,2); 
for x=1:365
    yearly_SImean_NYA(x,1)=x;
    yearly_SImean_NYA(x,2)=mean(NYA_SI(NYA_SI(:,1)==x,2),'omitnan');
    yearly_SIstd_NYA(x,1)=x;
    yearly_SIstd_NYA(x,2)=std(NYA_SI(NYA_SI(:,1)==x,2),'omitnan');
end

figure
plot(yearly_SImean_NYA(:,1),yearly_SImean_NYA(:,2),'k','LineWidth',2)
hold on
errorbar(yearly_SImean_NYA(:,1),yearly_SImean_NYA(:,2),yearly_SIstd_NYA(:,2),'k','CapSize',2);
xlim([0 365])
sumdoem=cumsum(eomday(2011,1:12));
xticks([0 sumdoem(1:end-1)]);
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel('Sea Ice Area Fraction / %')
title('Zeppelin Mountain Daily Mean \pm \sigma  (2011-2014)')
ax=gca;
ax.FontWeight='bold';
ax.FontSize=15;

yyaxis right
plot(yearly_HGmean_NYA(:,1),yearly_HGmean_NYA(:,2),'LineWidth',1.5)
errorbar(yearly_HGmean_NYA(:,1),yearly_HGmean_NYA(:,2),yearly_HGstd_NYA(:,2),'CapSize',2);
ylabel('Hg / ng m^{-3}')
hold off

%% Plotting daily SI median for VRS
% Formatting, averaging, and plotting data for daily mean for the period
SI_VRS.(1)=datenum(SI_VRS.(1));
VRS_SI=table2array(SI_VRS);
doy_vrs_si=date2doy(VRS_SI(:,1));
doy_vrs_si=floor(doy_vrs_si);
VRS_SI(:,1)=doy_nya_si;
 
yearly_SImean_VRS=NaN(365,2);
yearly_SIstd_VRS=NaN(365,2);
for x=1:365
    yearly_SImean_VRS(x,2)=mean(VRS_SI(VRS_SI(:,1)==x,2),'omitnan');
    yearly_SImean_VRS(x,1)=x;
    yearly_SIstd_VRS(x,1)=x;
    yearly_SIstd_VRS(x,2)=std(VRS_SI(VRS_SI(:,1)==x,2),'omitnan');
end

figure
plot(yearly_SImean_VRS(:,1),yearly_SImean_VRS(:,2),'k','LineWidth',2)
hold on
errorbar(yearly_SImean_VRS(:,1),yearly_SImean_VRS(:,2),yearly_SIstd_VRS(:,2),'k','CapSize',2);
xlim([0 365])
xticks([0 sumdoem(1:end-1)]);
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})

ylabel('Sea Ice Area Fraction / %')
title('Villum Research Station Daily Mean \pm \sigma  (2011-2014)')
ax=gca;
ax.FontWeight='bold';
ax.FontSize=15;

yyaxis right
plot(yearly_HGmean_VRS(:,1),yearly_HGmean_VRS(:,2),'LineWidth',1.5)
errorbar(yearly_HGmean_VRS(:,1),yearly_HGmean_VRS(:,2),yearly_HGstd_VRS(:,2),'CapSize',2);
ylabel('Hg / ng m^{-3}')
hold off

%% Sorting Data for AMDEs
% Indexing for the SI data for AMDEs
SI_VRS.(1)=datenum(SI_VRS.(1));
SI_VRS=table2array(SI_VRS);
SI_VRS(:,1)=floor(SI_VRS(:,1));
AMDE_VRS_11_14=table2array(AMDE_VRS_11_14);

AMDE_VRS_SI=NaN(length(AMDE_VRS_11_14),2);
AMDE_VRS_SI(:,1)=AMDE_VRS_11_14(:,1);
for x=1:length(AMDE_VRS_11_14)
    AMDE_VRS_SI(x,2)=SI_VRS(SI_VRS(:,1)==AMDE_VRS_11_14(x,1),2);
end

SI_NYA.(1)=datenum(SI_NYA.(1));
SI_NYA=table2array(SI_NYA);
SI_NYA(:,1)=floor(SI_NYA(:,1));

AMDE_NYA_SI=NaN(length(AMDE_NYA11_14),2);
AMDE_NYA_SI(:,1)=AMDE_NYA11_14(:,1);
for x=1:length(AMDE_NYA11_14)
    AMDE_NYA_SI(x,2)=SI_NYA(SI_NYA(:,1)==AMDE_NYA11_14(x,1),2);
end

%% Loading, Cleaning, Filtering, and Merging Met Data at VRS
MetVRS14=readtable('StNordTempRH2014_18.xlsx','ReadVariableNames',true);
MetVRS14 = removevars(MetVRS14, {'Station','x112','x113','x122','x123','x301','x305','x365','x371','x401','x504','x550','x601','x603','x609','x801'});
MetVRS14.time=datetime(MetVRS14.(1),MetVRS14.(2),MetVRS14.(3),MetVRS14.(4),0,0);
MetVRS14=table(MetVRS14.time,MetVRS14.(5),MetVRS14.(6));

MetVRS14=table2timetable(MetVRS14);
MetVRS14=retime(MetVRS14,'daily','mean');
MetVRS14=timetable2table(MetVRS14);

MetVRS14=MetVRS14(year(MetVRS14{:,1})==2014,:);
MetVRS14.(1)=datenum(MetVRS14.(1));
MetVRS14=table2array(MetVRS14);

AMDE_VRS_14=AMDE_VRS_11_14(AMDE_VRS_11_14(:,2)==2014,:);

VRS_Met_11_13=readtable('VRS_Met_11_14.xlsx','ReadVariableNames',false);

VRS_Met_14=NaN(length(AMDE_VRS_14),3);
VRS_Met_14(:,1)=AMDE_VRS_14(:,1);
for x=1:length(AMDE_VRS_14)
    VRS_Met_14(x,2:3)=MetVRS14(MetVRS14(:,1)==AMDE_VRS_14(x,1),2:3);
end

VRS_Met_14=array2table(VRS_Met_14);
VRS_Met_14.(1)=datetime(VRS_Met_14.(1),'ConvertFrom','datenum');
VRS_Met_11_13.Properties.VariableNames={'DateTime','Temp','RH'};
VRS_Met_14.Properties.VariableNames=VRS_Met_11_13.Properties.VariableNames;
VRS_Met_11_14=[VRS_Met_11_13;VRS_Met_14];

%% Loading, Merging, and Calcualting Stats for all parameters at VRS

addpath('C:\Users\au595124\Documents\PhD Thesis\Courses\Abisko')
VRS_Rad=readtable('AMDE_RAD_VRS.csv');
VRS_Rad=removevars(VRS_Rad,'Var1');

VRS_Met_11_14.Rad=VRS_Rad.rad;
VRS_Met_11_14.SI=AMDE_VRS_SI(:,2);
VRS_Met_11_14.Month=month(VRS_Met_11_14.(1));
statsVRS=grpstats(VRS_Met_11_14,{'Month','Temp','RH','Rad','SI'},[],'DataVars',{'Month','Temp','RH','Rad','SI'});
month_VRS_met_all=statsVRS{:,{'Month','Temp','RH','Rad','SI'}};

Monthly_avg_VRS(:,1)=[2:5];
for x=2:5
    for y=1:5
        Monthly_avg_VRS(x-1,y)=nanmean(month_VRS_met_all(month_VRS_met_all(:,1)==x,y));
    end
end

%% Loading, Merging, and Calcualting Stats for all parameters at NYA
NYA_Met_11_14=readtable('NYA_Met_11_14.xlsx','ReadVariableNames',false);
NYA_Met_11_14 = removevars(NYA_Met_11_14, {'Var2','Var4','Var5','Var6','Var7','Var8','Var11'});
NYA_Met_11_14.Properties.VariableNames={'Time','SWD','Temp','RH'};
NYA_Met_11_14.SI=AMDE_NYA_SI(:,2);
NYA_Met_11_14.Month=month(NYA_Met_11_14.Time);
statsNYA=grpstats(NYA_Met_11_14,{'Month','Temp','RH','SWD','SI'},[],'DataVars',{'Month','Temp','RH','SWD','SI'});
month_NYA_met_all=statsNYA{:,{'Month','Temp','RH','SWD','SI'}};

Monthly_mean_NYA(:,1)=[1:6];
for x=1:6
    for y=1:5
        Monthly_mean_NYA(x,y)=nanmean(month_NYA_met_all(month_NYA_met_all(:,1)==x,y));
    end
end
Monthly_mean_NYA(2,:)=[];
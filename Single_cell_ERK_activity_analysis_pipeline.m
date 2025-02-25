%% Single cell ERK activity analysis
% Julio Cesar Sanchez Rendon, PhD student, University of Tübingen
% 19/02/2025
clc, clear, close all

%% Track extraction and formating from TrackMate (xml) files
% Call tracks file
tracks_path = '/Users/Julio/Desktop/ERK analysis example';
tracks_file = 'Example_CFP.xml';
time_step = 10; % [min] time between frames
fcal = 0.28;   % [um/pixel] Microscopy calibration factor
img_size =  1608; % [pixel] Size of original image
Max_distance = 30; % [pixel] Maximum distance of frame-to-frame linking
tracks = readstruct([tracks_path,filesep,tracks_file]); 
nFrames = tracks.Settings.ImageData.nframesAttribute;
img_height = tracks.Settings.ImageData.heightAttribute;

% Clean short tracks
tracks = tracks.Model.AllTracks.Track;
Track_size = cell2mat({tracks.NUMBER_SPOTSAttribute});
Track_exclude = find(Track_size<0.8*nFrames);
tracks(Track_exclude) = [];

% Save individual track coordinates and time point per nucleus
track_nuclei = {};
time = [];
img_center = img_size/2;
j = 1;
for j=1:size(tracks,2)
    single_track = struct2table(tracks(j).Edge);
    single_track = single_track(:,7:9);
    single_track = table2array(sortrows(single_track,1));
    % Track info: time, x-coor, y-coor
    track_nuclei{1,j} = single_track;
    % correction for gap or duplication in time step
    track_nuclei{1,j}(~mod(track_nuclei{1,j}(:,1)+0.5,1)==0,:) = [];
    % Correction for tracks with 3 or more split events: trim coordinates
    % for duplicated time step
    if tracks(j).NUMBER_SPLITSAttribute > 1
        u_times = tabulate(track_nuclei{1,j}(:,1));
        str_track_end_point = find(track_nuclei{1,j}(:,1) == u_times(find(u_times(:,2)==1,1,'last'),1));
        track_nuclei{1,j}((str_track_end_point+1):end,:) = [];
        single_track = track_nuclei{1,j};
    end

    % Track properties
    track_nuclei{2,j} = sqrt(sum(diff((single_track(:,2:3)).*fcal).^2,2));% Displacement [um]
    track_nuclei{3,j} = track_nuclei{2,j}.*(60/time_step); % speed [um/h]
    track_nuclei{4,j} = sum(track_nuclei{2,j}); % Total traveled distance [um]
    track_nuclei{5,j} = (sqrt(sum((single_track(end,2:3)-single_track(1,2:3)).^2))*fcal)./track_nuclei{4,j}; % Persistence
    A = diff(single_track(:,2:3)); % x and y vector-coor for nuclei displacements
    B = single_track(2:end,2:3) - repmat(img_center,length(single_track(2:end,1)),2); % x and y coor w.r.t center
    track_nuclei{6,j} = atan2d((A(:,1).*B(:,2)-A(:,2).*B(:,1)),(A(:,1).*B(:,1)+A(:,2).*B(:,2))); % Angle w.r.t center of image [0 to +/-180°]
    Cen_coor_ang = atan2d(B(:,2),B(:,1));
    track_nuclei{7,j} = tracks(j).NUMBER_SPLITSAttribute;
    single_track = [];
    A = [];
    B = [];
    Cen_coor_ang = [];
    time = [time;track_nuclei{j}(:,1)];
end
time = unique(time);
time(mod(time-round(time),0.5)~=0) = [];
time((time-round(time)) ~= -0.5) = [];

% Correction for tracks with 2 or more split events
i = 1;
for i=find(table2array(table(track_nuclei{7,:}))>=2)
    c = find(ismember(track_nuclei{1,i}(:,1),setdiff(track_nuclei{1,i}(:,1),time)));
    track_nuclei{1,i}(c,:) = []; 
    u_times = tabulate(track_nuclei{1,i}(:,1));
    str_track_end_point = find(track_nuclei{1,i}(:,1) == u_times(find(u_times(:,2)==1,1,'last'),1));
    track_nuclei{1,i}((str_track_end_point+1):end,:) = []; 
end

%% ERK intensity and single cell track matching
% Call images
ERK_img = 'Example_ERK.tif';     % ERK intensity image name
seg_img = 'MASK_Original_segmented_nuclei.tif'; % Nuclei segmentation image name
ERK_img_num = imfinfo([tracks_path,filesep,ERK_img]);
ERK_img_num = length(ERK_img_num);

% Loop for time steps
i = 1;
for i=1:(ERK_img_num-1)
    % Call segmented nuclei mask
    bin_img = imread([tracks_path,filesep,seg_img],i);
    bin_img(bin_img == 255) = 1; % Threshold set for 8 bit image
    % Call ERK intensity image
    ERK_int_img = imread([tracks_path,filesep,ERK_img],i);
    ERK_int_img = ERK_int_img.*double(bin_img);
    % Calculate labeled matrix from segmented mask
    label_img = bwlabel(bin_img);
    % Match nuclei centroid from track with labeled matrix
    j = 1;
    for j = 1:length(track_nuclei)
        coor_xy_t = track_nuclei{1,j}(track_nuclei{1,j}(:,1) == double(time(i)),2:3);
        if ~isempty(coor_xy_t)
            n_xy = round(coor_xy_t,TieBreaker="tozero");
            n_xy(n_xy < 1) = 1;
            n_xy(n_xy > size(ERK_int_img,1)) = size(ERK_int_img,1);
            ERK_int_aux = ERK_int_img.*(label_img == label_img(n_xy(2),n_xy(1))); % x-coordinate is column & y-coordinate is row
            ERK_int_aux(ERK_int_aux == 0) = NaN;clc
            Nuclei_ERK_int(j,i) = median(ERK_int_aux,"all","omitnan");
        else
            coor_xy_t = [NaN,NaN];
            Nuclei_ERK_int(j,i) = NaN;
        end
    end
    Nuclei_ERK_int(isnan(Nuclei_ERK_int(:,i)),i) = 0;
end

%% Save results
% Store ERK intensity with nuclei tracks
i=1;
for i =1:length(track_nuclei(1,:))
nuc_erk = Nuclei_ERK_int(i,:)';
if length(nuc_erk) < length(track_nuclei{1,i}(:,1))
nuc_erk(end+1:length(track_nuclei{1,i}(:,1)),1) = NaN;
elseif length(nuc_erk) > length(track_nuclei{1,i}(:,1))
nuc_erk = nuc_erk(1:length(track_nuclei{1,i}(:,1)));
end
track_nuclei{1,i}(:,4) = nuc_erk;
end
vars = {'Track','Displacement','Speed','Total traveled distance','Persistance','Angle w.r.t center','Spliting events'}';
track_nuclei = cat(2,vars,track_nuclei);

% save analysis in cell per nucleus: time, x-coor, y-coor, ERK intensity
save([tracks_path,filesep,'Analysis.mat'],'track_nuclei','Nuclei_ERK_int','time_step','img_size','fcal','tracks_path')

%% Single cell ERK activity analysis figures 
% Heatmap of ERK activity organized by signal correlation 
figure(1)
htmp_type = "position";
switch htmp_type
    case "correlation"
        [H,~,Tree_order]=dendrogram(linkage(squareform(pdist(Nuclei_ERK_int)),"average",'correlation'),0,'Orientation','left'); % order by correlation
    
    case "position"
        Nuc_dist_to_center = [];
        for i = 1:length(track_nuclei(1,2:end))
            last_pos = length(track_nuclei{1,i+1}(:,4));
            int_nan = find(isnan(track_nuclei{1,i+1}(:,4)));
            int_zero = find(track_nuclei{1,i+1}(:,4) == 0);
            if (last_pos-(length(int_nan) + length(int_zero)))/last_pos < 0.85
                %last_pos = int_nan;
                continue
            end
            Nuc_dist_to_center(i,:) = [sqrt(sum((track_nuclei{1,i+1}(last_pos,2:3)-(img_size/2)*fcal).^2)),i];
            last_pos = [];
            int_nan = [];
            int_zero = [];
        end
        Tree_order = sortrows(Nuc_dist_to_center,1);
        Tree_order = Tree_order(:,2);
        Tree_order(Tree_order == 0) = [];
end
hmp = heatmap(Nuclei_ERK_int(Tree_order,:),'GridVisible','off');
time_plot = [1:size(Nuclei_ERK_int,2)]*time_step/60;
timeLabel = string(time_plot);
timeLabel(mod(time_plot,2)~=0) = '';
hmp.XLabel = 'Time (h)';
hmp.YLabel = 'Single cell tracks';
hmp.XDisplayLabels = timeLabel;
hmp.YDisplayLabels = repmat({''},1,length(hmp.YData));
hmp.Colormap = jet;
hmp.ColorLimits = [0.6 1.3];
title(['Total number of cell tracks = ',num2str(length(Tree_order))])
set(gca,'FontName','Arial','FontSize',16)
set(gcf,"Color",'w','Units','centimeters',"Position",[0,0,40,40])
print([tracks_path filesep 'heatmap.eps'],'-depsc')

% Single-cell ERK signal statistics 
figure(2)
ERK_time_series = Nuclei_ERK_int(Tree_order,:);
median_erk_time_series = median(ERK_time_series(:,1:length(time_plot)));
iqr_erk_time_series = quantile(ERK_time_series(:,1:length(time_plot)),[0.25 0.75]);
plot(time_plot,median_erk_time_series,'-b','LineWidth',1)
hold on
patch([time_plot,flip(time_plot)],[iqr_erk_time_series(1,:),flip(iqr_erk_time_series(2,:))],'b','FaceAlpha',0.1,'EdgeColor','w')
hold off
box off
xlim([time_plot(1),time_plot(end)])
ylim([0.4 1.8])
ylabel('Median of ERK signal')
xlabel('Time (h)')
set(gca,'FontName','Arial','FontSize',16)
set(gcf,"Color",'w','Units','centimeters',"Position",[0,0,40,25])
print([tracks_path filesep 'heatmap_statistics.pdf'],'-dpdf','-bestfit')

% Plot of radial aligment distribution
figure(3)
i=1;
Radial_align = [];
for i=2:length(track_nuclei(6,2:end))
    angles = track_nuclei{6,i}(1:end);
    Radial_align = [Radial_align;angles];
end
polarhistogram(deg2rad((abs(Radial_align))),24);
ph_ax = gca;
ph_ax.ThetaLim = [0 180];
set(gcf,'color','w')
print([tracks_path filesep 'Radial_alignment.pdf'],'-dpdf')

% Plot of cell speed distribution
figure(4)
cell_speeds = cell2mat(track_nuclei(3,2:end)');
cell_speeds(cell_speeds>Max_distance*fcal*60/time_step) = []; % Remove speeds with displacement higher than maximum linked distance
i=1;
for i = 1:3
    cell_speeds = rmoutliers(cell_speeds,"median"); % remove outliers iteratively
end
boxplot(cell_speeds,'Symbol','o');
ylim([0 35])
ylabel('Cell speed (µm/h)')
set(findobj(gca,'type','line'),'linew',2)
set(gca,'FontName','Arial','FontSize',16,'XTick',[])
set(gcf,'color','w')
print([tracks_path filesep 'Cell_speed.pdf'],'-dpdf')

% Plot of cell directionality distribution
figure(5)
boxplot(cell2mat(track_nuclei(5,2:end)),'Symbol','o')
ylabel('Directionality (a.u.)')
ylim([0 1])
set(findobj(gca,'type','line'),'linew',2)
set(gca,'FontName','Arial','FontSize',16,'XTick',[])
set(gcf,'color','w')
print([tracks_path filesep 'Cell_directionality.pdf'],'-dpdf')

% Generation of CVS file with radial alignment angles, cell speeds and linear persistance
Info_cell = {deg2rad((abs(Radial_align))),cell_speeds,cell2mat(track_nuclei(5,2:end))'};
length_cell = cellfun(@length,Info_cell);
max_size = max(length_cell);
i = 1;
for i = 1:3
    if length_cell(i) < max_size
        Info_cell{i} = [Info_cell{:,i};NaN(max_size-length_cell(i),1)];
    end
end

Migration_analysis = table(Info_cell{1},Info_cell{2},Info_cell{3},...
    'VariableNames',{'Radial_alignment','Cell_speed','Linear_persistence'});

writetable(Migration_analysis,[tracks_path,filesep,'Migration_analysis.csv'],'Delimiter','tab','WriteVariableNames',true);


%% ARCOS formating of tracks and ERK signal activation definition
ARCOS_matrix = [];
i = 1;
ERK_thrhld = 0.1; % default: 0.02 Signal threshold for ERK signal activation
Median_global = median(Nuclei_ERK_int);

for i=1:length(track_nuclei(1,2:end))
    AM_aux = track_nuclei{1,i+1};
    AM_aux(isnan(AM_aux(:,4)),:) = [];
    AM_aux(:,5) =  i;
    % remove duplicated time points
    [~,unique_tp] = unique(AM_aux(:,1));
    AM_aux = AM_aux(unique_tp,:);
    % remove points beyond 400 time steps
    AM_aux(AM_aux(:,1)+0.5>400,:) = [];
    
    % Filter by positive values after long-range median filter (30 steps)
    % Removal of global effect on indiviual signal
    AM_aux(:,4) = AM_aux(:,4) - Median_global(AM_aux(:,1)+0.5)';
    % Removal of trend in individual signal
    ERK_int_filter = smoothdata(AM_aux(:,4),'movmedian',30);
    AM_aux(:,4) = AM_aux(:,4) - ERK_int_filter;
    AM_aux(:,4) = AM_aux(:,4) > ERK_thrhld;

    ARCOS_matrix = [ARCOS_matrix;AM_aux];
    AM_aux = [];
end

% Formating of ERK activation matrix for ARCOS algorithm
ARCOS_matrix = sortrows(ARCOS_matrix,1);
ARCOS_matrix(:,1) = round(ARCOS_matrix(:,1),0,TieBreaker="tozero");
DBscan_dist = ARCOS_matrix(ARCOS_matrix(:,1)==0,2:3);

% Video of 50 frames with activated and unactivated single cell nuclei
length_video = min([50,size(Nuclei_ERK_int,2)]);
figure
v = VideoWriter([tracks_path,filesep,'ERK_waves_median'],'MPEG-4');
open(v)
for i = 1:length_video
    test_mat = ARCOS_matrix(ARCOS_matrix(:,1)==(i-1),:);
    scatter(test_mat(:,2),test_mat(:,3),30,test_mat(:,4),'filled'); colormap(viridis); set(gca,'YDir','reverse'); xlim([0 1608]); ylim([0 1608]); axis square
    text(50,50,['Frame:',' ',num2str(i-1)])
    frame = getframe(gcf);
    writeVideo(v,frame)
end
close(v)

% Saving of ERK activation matrix in csv format for easy reading in R 
ARCOS_matrix = array2table(ARCOS_matrix,'VariableNames',{'t','x','y','m','id'});
writetable(ARCOS_matrix,[tracks_path filesep 'ARCOS_matrix.csv'])

% Calculation of DBSCAN search distance from distance ditribution of nuclei
DBscan_dist = pdist2(DBscan_dist,DBscan_dist,'euclidean','smallest',5);
DBscan_dist = median(DBscan_dist(end,:));
fprintf('Search radius for DBSCAN algorithm: %g \n',DBscan_dist)

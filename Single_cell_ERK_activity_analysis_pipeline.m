%% Single cell ERK activity analysis
% Julio Cesar Sanchez Rendon, PhD student, University of Tübingen
% 09/08/2024
clc, clear, close all

%% Track extraction and formating from TrackMate (xml) files
% Call tracks file
tracks_path = '/Users/Julio/Desktop/20240801_LM_ERK_nuclei/Ctr_4';
tracks_file = 'Ctr_IF5_Pos3_CFP_16h.xml';
time_step = 5; % [min] time between frames
fcal = 0.28;   % [um/pixel] Microscopy calibration factor
img_size =  1608; % [pixel] Size of original image
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
    track_nuclei{2,j} = sqrt(sum(diff(single_track(:,2:3)).^2,2)).*fcal;% Displacement [um]
    track_nuclei{3,j} = track_nuclei{2,j}./(time_step/60); % speed [um/h]
    track_nuclei{4,j} = sum(track_nuclei{2,j}); % Total traveled distance [um]
    track_nuclei{5,j} = (sqrt(sum((single_track(end,2:3)-single_track(1,2:3)).^2))*fcal)./track_nuclei{4,j}; % Persistance
    A = diff(single_track(:,2:3)); % x and y vector-coor for nuclei displacements
    B = single_track(2:end,2:3) - repmat(img_center,length(single_track(2:end,1)),2); % x and y coor w.r.t center
    track_nuclei{6,j} = atan2d((A(:,1).*B(:,2)-A(:,2).*B(:,1)),(A(:,1).*B(:,1)+A(:,2).*B(:,2))); % Angle w.r.t center of image [0 to +/-180°]
    Cen_coor_ang = atan2d(B(:,2),B(:,1));
    track_nuclei{7,j} = (A(:,1).*cosd(Cen_coor_ang) + A(:,2).*sind(Cen_coor_ang))./(time_step/60); % Radial strain rate [um/h]: traveled radial distance/time_step
    track_nuclei{8,j} = tracks(j).NUMBER_SPLITSAttribute;
    single_track = [];
    A = [];
    B = [];
    Cen_coor_ang = [];
    time = [time;track_nuclei{j}(:,1)];
end
time = unique(time);
time(mod(time-round(time),0.5)~=0) = [];
time((time-round(time)) ~= -0.5) = [];

% Correction for tracks with 3 or more split events
i = 1;
for i=find(table2array(table(track_nuclei{8,:}))>=3)
    c = find(ismember(track_nuclei{1,i}(:,1),setdiff(track_nuclei{1,i}(:,1),time)));
    track_nuclei{1,i}(c,:) = []; 
    u_times = tabulate(track_nuclei{1,i}(:,1));
    str_track_end_point = find(track_nuclei{1,i}(:,1) == u_times(find(u_times(:,2)==1,1,'last'),1));
    track_nuclei{1,i}((str_track_end_point+1):end,:) = []; 
end

%% ERK intensity and single cell track matching
% Call images
ERK_img = 'Ctr_IF5_Pos3_ERK_16h.tif';     % ERK intensity image name
seg_img = 'MASK_Nuclei_segmentation.tif'; % Nuclei segmentation image name
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
vars = {'Track','Displacement','Speed','Total traveled distance','Persistance','Angle w.r.t center','Radial strain rate','Spliting events'}';
track_nuclei = cat(2,vars,track_nuclei);

% save analysis in cell per nucleus: time, x-coor, y-coor, ERK intensity
save([tracks_path,filesep,'Analysis.mat'],'track_nuclei','Nuclei_ERK_int','time_step','img_size','fcal','tracks_path')

%% ARCOS formating of tracks and ERK signal activation definition
ARCOS_matrix = [];
i = 1;
ERK_thrhld = 0.01; % default: 0.02 Signal threshold for ERK signal activation
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
v = VideoWriter('ERK_waves_median','MPEG-4');
open(v)
for i = 1:50
    test_mat = ARCOS_matrix(ARCOS_matrix(:,1)==i,:);
    scatter(test_mat(:,2),test_mat(:,3),30,test_mat(:,4),'filled'); colormap(viridis); set(gca,'YDir','reverse'); xlim([0 1608]); ylim([0 1608]); axis square
    text(50,50,['Frame:',' ',num2str(i)])
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
fprintf('Search radious for DBSCAN algorithm: %g \n',DBscan_dist)


%% Figures 
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

% Plot of radial aligment
figure(2)
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
print([tracks_path filesep 'Radial_aligment.pdf'],'-dpdf')

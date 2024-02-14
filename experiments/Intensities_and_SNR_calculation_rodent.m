% Analysis script for intensities and SNR calculation of Rat Doppler data

% Rabut, Norman, Griggs et al., Science Translational Medicine (2024)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%           Locate Data             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Doing this for a single representative rodent
% Session - 20221008

current_dir = cd;

Conditions = cell(1, 5);
for i = 1:5
    Conditions{1,i} = w2b.io.getRodentFilePath('project_record_index', 20+i, ...
        'suppress_screen_output', true);
end

names = {'No implant';'1mm';'2mm';'3mm';'Mesh'};


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      Display Power Doppler         %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
for cond = 1:length(Conditions)
    clear Dop
    path = [Conditions{1,cond}];
    cd(path)
    load('Dop.mat')
    subplot(2,3,cond)
    imagesc(sqrt(squeeze(mean(Dop,3))))
    colormap hot
    title(names(cond))
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Resize Doppler images to all be at same depth   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Doing this for a single representative rodent
% Session - 20221008
% Hard-coding specific depths to standardize depth

cond = 1 ;
clear Dop
path = [Conditions{1,cond}];
cd(path)
load('Dop.mat')
Conditions_reshape{1,cond}= Dop(29:end,:,:);

cond = 2 ;
clear Dop
path = [Conditions{1,cond}];
cd(path)
load('Dop.mat')
Conditions_reshape{1,cond}= Dop(27:end,:,:);

cond = 3 ;
clear Dop
path = [Conditions{1,cond}];
cd(path)
load('Dop.mat')
Conditions_reshape{1,cond}= Dop(27:end,:,:);

cond = 4 ;
clear Dop
path = [Conditions{1,cond}];
cd(path)
load('Dop.mat')
Conditions_reshape{1,cond}= Dop(25:end,:,:);

cond = 5 ;
clear Dop
path = [Conditions{1,cond}];
cd(path)
load('Dop.mat')
Conditions_reshape{1,cond}= Dop(31:end,:,:);

figure
for cond = 1:length(Conditions_reshape)
    subplot(2,3,cond)
    imagesc(sqrt(squeeze(mean(Conditions_reshape{1,cond},3))))
    colormap hot
    title(names(cond))
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      Compare total intensity at baseline state      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear intens
for cond = 1:length(Conditions_reshape)
    Dop = Conditions_reshape{1,cond};
    intens(cond) = squeeze(mean(mean(mean(Dop(:,:,1:30),1),2),3)) ; % Baseline = first 30 blocks
end
intens = intens/intens(1) ;

figure
plot(intens,'s','MarkerSize',10,'MarkerEdgeColor',[0.98,0.29,0.29])
title('Normalized Doppler intensity')
set(gca,'xtick',[1:length(intens)],'xticklabel',names(1:length(intens)))



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    Compare SNR at two depths (at baseline state)    %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear SNR_cort SNR_deep_struct CNR_cort
for cond = 1:length(Conditions_reshape)
    Dop = Conditions_reshape{1,cond};

    clear cortex_region deep_struct_region
    cortex_region = squeeze(Dop(4:8,48:58,2));
    deep_struct_region = squeeze(Dop(32:35,51:77,2));

    i=1;
    for d = 1:size(cortex_region,1)
        cortex_line2study = squeeze(cortex_region(d,:));
        loc_max_cortex = islocalmax(cortex_line2study);
        loc_min_cortex = islocalmin(cortex_line2study);
        x = [1:length(cortex_line2study)] ;
        SNR_cort(i,cond) = mean(cortex_line2study(loc_max_cortex))/mean(cortex_line2study(loc_min_cortex));
        i = i+1 ;
        clf
    end

    i=1;
    for d = 1:size(deep_struct_region,1)
        deep_struct_line2study = squeeze(deep_struct_region(d,:));
        loc_max_deep_struct = islocalmax(deep_struct_line2study);
        loc_min_deep_struct = islocalmin(deep_struct_line2study);
        x = [1:length(deep_struct_line2study)] ;
        SNR_deep_struct(i,cond) = mean(deep_struct_line2study(loc_max_deep_struct))/mean(deep_struct_line2study(loc_min_deep_struct)) ;
        i = i+1 ;
        clf
    end
end

% Plot SNR in function of implant scenario
SNR_cort_norm = SNR_cort/squeeze(mean(SNR_cort(:,1)));
SNR_cort_mean = squeeze(mean(SNR_cort_norm,1));
SNR_cort_std = squeeze(std(SNR_cort_norm,1));

SNR_deep_struct_norm = SNR_deep_struct/squeeze(mean(SNR_deep_struct(:,1)));
SNR_deep_struct_mean = squeeze(mean(SNR_deep_struct_norm,1));
SNR_deep_struct_std = squeeze(std(SNR_deep_struct_norm,1));

figure,
subplot(2,1,1)
errorbar(SNR_cort_mean,SNR_cort_std,'s','MarkerSize',10,'MarkerEdgeColor',[0.98,0.29,0.29])
title('SNR around cortex')
set(gca,'xtick',[1:size(SNR_cort,2)],'xticklabel',names(1:size(SNR_cort,2)) )
axis([0.3 size(SNR_cort,2)+0.6 0 1.5])
subplot(2,1,2)
errorbar(SNR_deep_struct_mean,SNR_deep_struct_std,'s','MarkerSize',10,'MarkerEdgeColor',[0.55,0.07,0.07])
title('SNR deep in brain')
set(gca,'xtick',[1:size(SNR_cort,2)],'xticklabel',names(1:size(SNR_cort,2)))
axis([0.3 size(SNR_cort,2)+0.6 0 1.5])


% Plot SNR attenuation (in dB) in function of implant scenario
SNR_cort_db = mag2db(SNR_cort_mean);
SNR_cort_std_db = mag2db(SNR_cort_std);

SNR_deep_struct_db = mag2db(SNR_deep_struct_mean);
SNR_deep_struct_std_db = mag2db(SNR_deep_struct_std);

figure,
subplot(2,1,1)
errorbar(SNR_cort_db,SNR_cort_std_db,'s','MarkerSize',10,'MarkerEdgeColor',[0.98,0.29,0.29])
title('Attenuation in dB around cortex')
set(gca,'xtick',[1:size(SNR_cort,2)],'xticklabel',names(1:size(SNR_cort,2)) )
axis([0.3 size(SNR_cort,2)+0.6 -75 60])
subplot(2,1,2)
errorbar(SNR_deep_struct_db,SNR_deep_struct_std_db,'s','MarkerSize',10,'MarkerEdgeColor',[0.55,0.07,0.07])
title('Attenuation in dB deep in brain')
set(gca,'xtick',[1:size(SNR_cort,2)],'xticklabel',names(1:size(SNR_cort,2)))
axis([0.3 size(SNR_cort,2)+0.6 -75 60])


cd(current_dir);

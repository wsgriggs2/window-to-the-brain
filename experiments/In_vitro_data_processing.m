% Analysis script for in vitro Doppler

% Rabut, Norman, Griggs et al., Science Translational Medicine (2024)

%% Load doppler data previously processed with In_vitro_Power_Doppler_STM.mat script

data_path = w2b.io.getUserDataPath;

path_name = fullfile(data_path, 'in vitro', 'In_vitro_Doppler_Data_STM.mat');

load(path_name);

%% Creates cells with different conditions

names = {'No implant';'1 mm';'2 mm';'3 mm';'Mesh'} ;
Dopplers(1,:,:,:) = Dop_No ;
Dopplers(2,:,:,:) = Dop_1mm ;
Dopplers(3,:,:,:) = Dop_2mm ;
Dopplers(4,:,:,:) = Dop_3mm ;
Dopplers(5,:,:,:) = Dop_mesh ;

% Display Power Doppler images of all five conditions

figure
for implant = 1:5
    subplot(1,5,implant)
    imagesc(squeeze(sqrt(Dopplers(implant,30:250,:,15))))
    colormap hot
    title(names{implant})
end

%% Total intensities
intens = squeeze(mean(mean(Dopplers,2),3))./squeeze(mean(mean(mean(Dopplers(1,:,:,:),2),3),4));

figure
for rep = 1:15
    plot(intens(:,rep),'LineStyle','none','Marker',"o",'MarkerSize',12,'MarkerEdgeColor','none','MarkerFaceColor','k')
    hold on
end
axis([0.5 5.5 0.8 1.1])
title('Power Doppler intensity')
categories = {'No implant', '1 mm PMMA', '2 mm PMMA', '3 mm PMMA', 'Ti Mesh'};
xticks(1:5)
xticklabels(categories)

%% Calculate SNR in each column and at each depth
clear column depth

column{1} = [1:35] ;
column{2} = [40:75] ;
column{3} = [75:128] ;

depth{1} = [62:73] ;
depth{2} = [110:125] ;
depth{3} = [162:175] ;
depth{4} = [212:230] ;

for implant = 1:length(names)            % loop over all the implant conditions
    for c = 1:length(column)             % loop over the three columns of doppler tubings
        for d = 1:length(depth)          % loop over the four lines of doppler tubings (in depth)
            for rep = 1:size(Dopplers,4) % loop over the 15 replicates

                % Find the tubing
                clear Dop2study line2study gauss MaxAmpl signal noise
                Dop2study = squeeze(Dopplers(implant,depth{d}, column{c},rep)) ;
                [z,x] = find(Dop2study == (max(Dop2study(:))));
                if z < 4 ,Dop2study = squeeze(Dopplers(implant,depth{d}(1)-9:depth{d}(end), column{c},rep)) ; end
                [z,x] = find(Dop2study == (max(Dop2study(:))));
                
                ma = 1;
                for z_search = z-1 : z+1
                    line2study = Dop2study(z_search,:) ;
                    gauss = fit([1:length(line2study)]',line2study','gauss1');
                    
                    MaxAmpl(ma) = gauss.a1 ;
                    ma = ma + 1 ;
                end
                
                % Calculate SNR
                
                % signal
                signal = max(MaxAmpl(:)) ;
                
                % noise
                clear z_max x_max noise noise_window
                z_max = z - 2 + find(MaxAmpl == (max(MaxAmpl(:))));
                [zx, x_max] = find(Dop2study == (max(Dop2study(z_max,:))) );
                
                for x = 3:size(Dop2study,2)-2
                    noise_window(x-2) = squeeze(mean(Dop2study(z_max,x-2:x+2))) ; % sliding window of 5 pixels
                end
                noise = min(noise_window) ; % find the darkest window at this depth
                x_noise =  find(noise_window == noise ) +2;
                
                SNR(implant,d,c,rep) = signal/noise ;
                
            end 
        end
    end
end

%% Plot SNR in function of depth and conditions (dB)

stats_SNR = reshape(SNR,[5 4 3*size(Dopplers,4)]) ;
mean_SNR = squeeze(mean(stats_SNR,3)) ;
SNR_dB = mag2db( bsxfun(@rdivide,stats_SNR,mean_SNR(1,1)));
colors = [0 0.4470 0.7410 ; 0.8500 0.3250 0.0980 ; 0.9290 0.6940 0.1250 ; 0.4940 0.1840 0.5560 ; 0.4660 0.6740 0.1880 ] ;
axis_implant = [0.8 1.8 2.8 3.8 ; 0.9 1.9 2.9 3.9 ; 1 2 3 4 ; 1.1 2.1 3.1 4.1 ; 1.2 2.2 3.2 4.2] ;


figure
for implant = 1:5
    for rep = 1:size(stats_SNR,3)
        plot(axis_implant(implant,:),SNR_dB(implant,:,rep),'LineStyle','none','Marker',"o",'MarkerSize',16,'MarkerEdgeColor','none','MarkerFaceColor',[colors(implant,:)])
        hold on
    end
end
axis([0.5 4.5 -30 20])
categories = {'14 mm', '24 mm', '34 mm', '14 mm', '44 mm'};
xticks(1:5)
xticklabels(categories)
xlabel('depth')
ylabel('Attenuation (dB)')
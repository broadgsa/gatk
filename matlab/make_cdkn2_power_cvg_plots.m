clear all
close all

need_to_create_metrics_file = false;

if ( need_to_create_metrics_file )
    production_names = {'FHS_P10_766_CDKN2'; ...
                           'FHS_P113_140_CDKN2'; ...
                           'FHS_P11_767_CDKN2'; ...
                           'FHS_P13_769_CDKN2'; ...
                           'FHS_P141_168_CDKN2'; ...
                           'FHS_P14_770_CDKN2'; ...
                           'FHS_P15_771_CDKN2'; ...
                           'FHS_P169_196_CDKN2'; ...
                           'FHS_P17_773_CDKN2'; ...
                           'FHS_P197_224_CDKN2'; ...
                           'FHS_P1_757_CDKN2'; ...
                           'FHS_P225_252_CDKN2'; ...
                           'FHS_P26_782_CDKN2'; ...
                           'FHS_P27_783_CDKN2'; ...
                           'FHS_P28_784_CDKN2'; ...
                           'FHS_P2_758_CDKN2'; ...
                           'FHS_P309_336_CDKN2'; ...
                           'FHS_P365_392_CDKN2'; ...
                           'FHS_P393_420_CDKN2'; ...
                           'FHS_P3_759_CDKN2'; ...
                           'FHS_P421_448_CDKN2'; ...
                           'FHS_P449_476_CDKN2'; ...
                           'FHS_P477_504_CDKN2'; ...
                           'FHS_P4_760_CDKN2'; ...
                           'FHS_P505_532_CDKN2'; ...
                           'FHS_P533_560_CDKN2'; ...
                           'FHS_P561_588_CDKN2'; ...
                           'FHS_P57_84_CDKN2'; ...
                           'FHS_P5_761_CDKN2'; ...
                           'FHS_P617_644_CDKN2'; ...
                           'FHS_P645_672_CDKN2'; ...
                           'FHS_P673_700_CDKN2'; ...
                           'FHS_P6_762_CDKN2'; ...
                           'FHS_P701_728_CDKN2'; ...
                           'FHS_P729_756_CDKN2'; ...
                           'FHS_P757_784_CDKN2'; ...
                           'FHS_P85_112_CDKN2'; ...
                           'FHS_P8_764_CDKN2'; ...
                          'FHS_P9_765_CDKN2';};
                      
    production_dir = '/humgen/gsa-hphome1/projects/FHS/production/analysis/';
    
    
    pilot_names = {'CEPH1_CDKN2'; 'CEPH2_CDKN2'; 'CEPH3_CDKN2'};
    pilot_dir = '/humgen/gsa-hphome1/projects/FHS/pilot/analysis/';
    
    
    power_ext = '_power.txt';
    cvg_ext = '_coverage.txt';
    
    
    production_cvg = [];
    production_power = [];
    pilot_cvg = [];
    
    % load production data
    
    for ii = 1:length(production_names)
        pow_file = strcat(production_dir,production_names{ii},power_ext);
        [chrompos, data, u1, u2, u3, u4] = textread(pow_file,'%s\t%f\t%f\t%f\t%f\t%f','headerlines',1);
        production_power = [production_power data(:,1)];
        cvg_file = strcat(production_dir,production_names{ii},cvg_ext);
        [chrompos,data] = textread(cvg_file,'%s\t%d','headerlines',1);
        production_cvg = [production_cvg data(:,1)];
    end
    
    % load pilot data
    
    for ii = 1:length(pilot_names)
        cvg_file = strcat(pilot_dir,pilot_names{ii},cvg_ext);
        [chrompos,data] = textread(cvg_file,'%s\t%d','headerlines',1);
        pilot_cvg = [pilot_cvg data(:,1)];
    end
    
    % grab the raw positions on chromosome 9
    
    pos = [];
    
    for ii = 1:length(chrompos)
        g = chrompos{ii};
        h = find(g == ':');
        pos = [pos; str2num(g(h+1:end))];
    end
    
    % import the amplicon list
    
    [amp_start, amp_end] = textread('/humgen/gsa-hphome1/projects/FHS/interval_lists/chr9_amplicons.interval_list','%d\t%d');
    
    % now make a plottable target_region variable
    % and an amplicon_start variable in that same space
    target = 0;
    target_region = zeros(size(pos));
    amplicon_start = zeros(size(pos));
    prevPos = pos(1) - 1;
    
    for ii = 1:length(pos)
        if ( pos(ii) == prevPos - 1 )
            target_region(ii) = target;
            target = target + 1;
        else
            target_region(ii) = target + 50;
            target = target + 51;
        end
        
        if ( any(pos(ii) == amp_start) )
            amplicon_start(ii) = 1; % yes we want this to be boolean
        end
        
    end
    
    
    % now save the data
    
    save('framingham_gene_CDKN2_metrics.mat','pilot_cvg','production_cvg','production_power','amplicon_start','target_region')
    
    
else


    load framingham_gene_CDKN2_metrics.mat
    
    % calculate median and mean
    mean_power = mean(production_power,2);
    mean_production_cvg = mean(production_cvg,2);
    mean_pilot_cvg = mean(pilot_cvg,2);
    median_power = median(production_power,2);
    median_production_cvg = median(production_cvg,2);
    median_pilot_cvg = median(pilot_cvg,2);
    
    % compute quartiles
    
    power_quartile1 = zeros(size(production_power(:,1)));
    power_quartile3 = zeros(size(power_quartile1));
    depth_prod_quartile1 = zeros(size(power_quartile1));
    depth_pilot_quartile1 = zeros(size(power_quartile1));
    depth_prod_quartile3 = zeros(size(power_quartile1));
    depth_pilot_quartile3 = zeros(size(power_quartile1));

    for ii = 1 : length(power_quartile1)
        power_quartile1(ii) = median(production_power(ii,find(production_power(ii,:) < median_power(ii)))')';
        power_quartile3(ii) = median(production_power(ii,find(production_power(ii,:) > median_power(ii)))')';
        depth_prod_quartile1(ii) = median(production_cvg(ii,find(production_cvg(ii,:) < median_production_cvg(ii)))')';
        depth_prod_quartile3(ii) = median(production_cvg(ii,find(production_cvg(ii,:) > median_production_cvg(ii)))')';
        depth_pilot_quartile1(ii) = median(pilot_cvg(ii,find(pilot_cvg(ii,:) < median_pilot_cvg(ii)))')';
        depth_pilot_quartile3(ii) = median(pilot_cvg(ii,find(pilot_cvg(ii,:) > median_pilot_cvg(ii)))')';
    end
    
    % take things into log space
    
    log_mean_power = log(1+mean_power);
    log_mean_production_cvg = log(1+mean_production_cvg);
    log_mean_pilot_cvg = log(1+mean_pilot_cvg);
    log_median_power = log(1+median_power);
    log_median_production_cvg = log(1+median_production_cvg);
    log_median_pilot_cvg = log(1+median_pilot_cvg);
    log_q1_power = log(1+power_quartile1);
    log_q3_power = log(1+power_quartile3);
    log_q1_production_cvg = log(1+depth_prod_quartile1);
    log_q3_production_cvg = log(1+depth_prod_quartile3);
    log_q1_pilot_cvg = log(1+depth_pilot_quartile1);
    log_q3_pilot_cvg = log(1+depth_pilot_quartile3);
    
    % get amplicon start positions
    
    amp_start = target_region(find(amplicon_start==1));

    % make plots
    
    grey = [0.7,0.7,0.7];
    
    h0 = figure;
    plot(target_region,mean_power,'r'), hold on
    plot(target_region,median_power,'k'), hold on
    plot(target_region,power_quartile1, 'color', grey), hold on
    plot(target_region, power_quartile3, 'color', grey), hold off
    set(gca,'xtick', amp_start)
    title('Power - Production')
    
    h1 = figure;
    plot(target_region, mean_pilot_cvg, 'r'), hold on
    plot(target_region, median_pilot_cvg, 'k'), hold on
    plot(target_region, depth_pilot_quartile1, 'color', grey), hold on
    plot(target_region, depth_pilot_quartile3, 'color', grey), hold off
    set(gca,'xtick', amp_start)
    title('Coverage - Pilot')
    
    h2 = figure;
    plot(target_region, mean_production_cvg, 'r'), hold on
    plot(target_region, median_production_cvg, 'k'), hold on
    plot(target_region, depth_prod_quartile1, 'color', grey), hold on
    plot(target_region, depth_prod_quartile3, 'color', grey), hold off
    set(gca,'xtick', amp_start)
    title('Coverage - Production')
    
    h3 = figure;
    plot(target_region,log_mean_power,'r'), hold on
    plot(target_region,log_median_power,'k'), hold on
    plot(target_region,log_q1_power, 'color', grey), hold on
    plot(target_region, log_q3_power, 'color', grey), hold off
    set(gca,'xtick', amp_start)
    title('Log power - Production')
    
    h4 = figure;
    plot(target_region, log_mean_pilot_cvg, 'r'), hold on
    plot(target_region, log_median_pilot_cvg, 'k'), hold on
    plot(target_region, log_q1_pilot_cvg, 'color', grey), hold on
    plot(target_region, log_q3_pilot_cvg, 'color', grey), hold off
    set(gca,'xtick', amp_start)
    title('Log coverage - Pilot')
    
    h5 = figure;
    plot(target_region, log_mean_production_cvg, 'r'), hold on
    plot(target_region, log_median_production_cvg, 'k'), hold on
    plot(target_region, log_q1_production_cvg, 'color', grey), hold on
    plot(target_region, log_q3_production_cvg, 'color', grey), hold off
    set(gca,'xtick', amp_start)
    title('Log coverage - Production')
    
end
close all
%load mean_coverage_per_exon_per_pool.mat <-- doesn't work from
%/sting/matlab -- but this does

imp = importdata('/humgen/gsa-hphome1/projects/FHS/results/production_averageDoCByExonAndGene.txt');
data = imp.data;
textdata = imp.textdata;

pilot_pos = ismember(textdata(:,2),'CEPH3')+ismember(textdata(:,2),'CEPH2')+ismember(textdata(:,2),'CEPH1');
pilot_data = data( find(pilot_pos == 1), : );
prod_data = data( find(pilot_pos == 0) , :);
pilot_text = textdata( find ( pilot_pos == 1) , : );
prod_text = textdata( find ( pilot_pos == 0 ), : );

lim = length(data);

% get unique names

gnames = [textdata(1,3)];
prevn = 'hello doofus';
for ii = 2 : lim/23  % first 1/23 of the file contains all the gene names (rest are repeats)
   n = textdata{ii,3};
   if ( length(prevn) == length(n) && prevn(1)==n(1) && prevn(end)==n(end) )
      % quick test to check if gene is same as previous
   else
       if ( ~ ismember(gnames,n) )
            gnames = [gnames; textdata(ii,3)];
            % more exhaustive test to see if gene name is novel
       end
       prevn = n;
   end
end

% plot the groups
num_genes = size(gnames)  % yes we want to print this

to_plot = 40;
plotno = 0;
filenamebaseprod = 'FHS_production_gene_cvg_boxplot_';
filenamebasepilot = 'FHS_pilot_gene_cvg_boxplot_';
for genes = 1 : to_plot : num_genes
    plotno = plotno + 1;
    pilot_positions = [];
    prod_positions = [];
    for g =  genes:genes+to_plot
        %n = gnames{g,1}
        if ( g < length(gnames) )
            pilot_positions = [pilot_positions; find(ismember(pilot_text(:,3),gnames{g,1}) == 1)];
            prod_positions = [prod_positions; find(ismember(prod_text(:,3),gnames{g,1}) == 1 )];
        end
    end
    depths_prod = prod_data(prod_positions,2);
    depths_pilot = pilot_data(pilot_positions,2);
    grenes_prod = prod_text(prod_positions,3);
    grenes_pilot = pilot_text(pilot_positions,3);
    h = figure
    hp = subplot(1,1,1)
    boxplot(depths_prod,grenes_prod,'plotstyle','compact','orientation','vertical','datalim',[0,10000],'whisker',0,'symbol','r+'), title('Production')
    set(hp,'YLim',[-1000,11000])
    y = figure
    yp = subplot(1,1,1)
    boxplot(depths_pilot,grenes_pilot,'plotstyle','compact','orientation','vertical','datalim',[0,10000],'whisker',0,'symbol','r+'), title('Pilot')
    set(yp,'YLim',[-1000,11000])
    
    % -- uncomment these lines to save the files -- 
    %saveas(h,strcat(filenamebaseprod,num2str(plotno)),'psc2');
    %saveas(y,strcat(filenamebasepilot,num2str(plotno)),'psc2');
end
        

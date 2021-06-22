clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colloid_folder = 'data/colloid/PsiW_neg4_PsiD_neg3/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colloid_cell_unload = matfile(join([colloid_folder,'colloid_cell.mat'],''));
colloid_cell = colloid_cell_unload.colloid_cell;
XYTcolloid_cell_unload = matfile(join([colloid_folder,'XYTcolloid_cell.mat'],''));
XYTcolloid_cell = XYTcolloid_cell_unload.XYTcolloid_cell;
colloid_log = matfile(join([colloid_folder,'log.mat'],''));
solute_log = matfile(join([colloid_log.solute_cell_folder,'log.mat'],''));
solute_cell_unload = matfile(join([colloid_log.solute_cell_folder,'solute_cell.mat'],''));
solute_cell = solute_cell_unload.solute_cell;
XYTsolute_cell_unload = matfile(join([colloid_log.solute_cell_folder,'XYTsolute_cell.mat'],''));
XYTsolute_cell = XYTsolute_cell_unload.XYTsolute_cell;



tic;

n_max = 0;
c_max = 0;


startyf = [];
startyf(1)=-1/solute_log.L_h;
for yValIndex = 1:length(XYTsolute_cell{1,2})
    if rem(yValIndex,2)==0
        startyf = [startyf XYTsolute_cell{1,2}(yValIndex)]; %#ok<AGROW>
    end
end
%startyf = XYTsolute_cell{1,2};
startxf = zeros(size(startyf));

stds = [];
counts = [];
% march in time
for count=1:length(XYTcolloid_cell{1,3})

    s = std(mean(colloid_cell{count,1},2));
    stds(count) = s;
    counts(count) = count;
    
end
toc;

plot(counts/0.481,stds,'-k','linewidth',1.5)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('time')
ylabel('std(n)')


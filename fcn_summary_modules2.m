function fcn_summary_modules2(ci,t1out_grp1,cmask,outprefix,t1out_grp2,summaryonly)
%%  --- Print and Compare Module Averaged Data --- %%
% INPUTS:
%   ci -            A module partition of the network nodes
%   t1out_grp1 -    Structure variable containing the output of tier1 for a group.
%   t1out_grp2 -    Structure variable containing the output of tier1 for group2.
%   cmask -         Binary matrix equal in size to data elements in the cellData
%                   that was used as tier1 input. It is used to isolate subsets of groups
%                   to be analyzed.
%   outprefix -     a prefix that is appended to all written out files
%
% Author: Evgeny Jenya Chumin (2023), Indiana University

%%
narginchk(2,6)
% check comm structure
 if ~isstruct(t1out_grp1)
    fprintf(2,'Tier1_out is not a structure. Exiting...')
    return
 end
% check node size vs size of community structure vector
if t1out_grp1.N ~= length(ci)
    fprintf(2,'Number of nodes does not match between inputs. Exiting...')
    return
end
% make sure communities are linearly indexed
if sum(ci==1)==0 || length(unique(diff(unique(ci)))) > 1 
    % if index doesnt start with one or if there are gaps in index
    ci = fcn_relabel_partitions(ci);
end
% optional cmask should be the same size as number of DATA cells in the
% original input cell structure
if exist('cmask','var')
    [r,c]=size(cmask);
    if r ~= t1out_grp1.Nr-1 || c ~= t1out_grp1.Nc-1
        fprintf(2,'Size of cmask does not match size data cells. Exiting...')
        return
    end
else
    cmask=ones(t1out_grp1.Nr-1,t1out_grp1.Nc-1);
end
if exist('t1out_grp2','var')
    if t1out_grp1.N ~= t1out_grp2.N || t1out_grp1.Nr ~= t1out_grp2.Nr || t1out_grp1.Nc ~= t1out_grp2.Nc
        fprintf(2,'Group1 and Group2 inputs do not match in size. Exiting...')
        return
    end
    grpF=2; % set group factor to two levels
else
    grpF=1;
end
if ~exist('summaryonly','var')
    summaryonly = 0;
end

%% Find cells from which to compute summary values
% pad zeros for row and column labels found in cell data
cmask = vertcat(zeros(1,t1out_grp1.Nc), horzcat(zeros(t1out_grp1.Nr-1,1),cmask));
% get row,col indices of included groups
[rw,cl]=find(cmask);
% decide on number of levels for anovan
rwF= length(unique(rw)); 
clF= length(unique(cl));
% find unique community labels
ciu = unique(ci);

% find largest group size
slength = max(max(t1out_grp1.S));
% text file output
gmeans = table;
if grpF==2
    gmeans_grp2 = table;
end
% initialize variables for anovan
alldata = cell(length(ciu),1);
for i=1:length(ciu)
    alldata{i,1}=double.empty;
end
if grpF==2
    allgrp = double.empty;
end
if rwF > 1
    allrw = double.empty;
end
if clF > 1
    allcl = double.empty;
end

% for every indexed group
for i=1:length(rw)
    d = t1out_grp1.cellData{rw(i),cl(i)};
    % for every community label (module)
    for cj=1:length(ciu)
        % extract the average and stdev regional value in community
        if sum(ci==ciu(cj))>1
            cmean = mean(d(ci==ciu(cj),:))';
            cstd = std(d(ci==ciu(cj),:))';
        else
            cmean = d(ci==ciu(cj),:)';
            cstd = nan(slength,1);
        end
        if summaryonly==0
       % concat data for anovan
        alldata{cj,1} = vertcat(alldata{cj,1},cmean);
        % update factors
        if cj == 1
            if grpF==2
                allgrp = vertcat(allgrp,ones(length(cmean),1));
            end
            if rwF > 1
                allrw = vertcat(allrw,zeros(length(cmean),1)+(rw(i)-1));
            end
            if clF > 1
                allcl = vertcat(allcl,zeros(length(cmean),1)+(cl(i)-1));
            end
        end
        end
        % label
        nm=[t1out_grp1.cellData{rw(i),1} '_' t1out_grp1.cellData{1,cl(i)} '_module' num2str(ciu(cj))];
        % pad with nan so it can be written into table
        if length(cmean)~=slength
            cmean(end+1:slength)=nan;
            cstd(end+1:slength)=nan;
        end
        % place data into table
        nm1=['mn_' nm];
        gmeans.(nm1)=cmean;
        nm2=['std_' nm];
        gmeans.(nm2)=cstd;
        clear cmean nm nm1 nm2
    end
    clear d Md
end
%% Run for group2 if one is given
if grpF==2
    % for every indexed group
    for i=1:length(rw)
        d = t1out_grp2.cellData{rw(i),cl(i)};
        % for every community label (module)
        for cj=1:length(ciu)
            % extract the average and stdev regional value in community
            if sum(ci==ciu(cj))>1
                cmean = mean(d(ci==ciu(cj),:))';
                cstd = std(d(ci==ciu(cj),:))';
            else
                cmean = d(ci==ciu(cj),:)';
                cstd = nan(slength,1);
            end
            if summaryonly==0
          % concat data for anovan
            alldata{cj,1} = vertcat(alldata{cj,1},cmean);
            % update factors
            if cj == 1
                if grpF==2
                    allgrp = vertcat(allgrp,ones(length(cmean),1)+1);
                end
                if rwF > 1
                    allrw = vertcat(allrw,zeros(length(cmean),1)+(rw(i)-1));
                end
                if clF > 1
                    allcl = vertcat(allcl,zeros(length(cmean),1)+(cl(i)-1));
                end
            end
            end
            % label
            nm=[t1out_grp2.cellData{rw(i),1} '_' t1out_grp2.cellData{1,cl(i)} '_module' num2str(ciu(cj))];
            % pad with nan so it can be written into table
            if length(cmean)~=slength
                cmean(end+1:slength)=nan;
                cstd(end+1:slength)=nan;
            end
            % place data into table
            nm1=['mn_' nm];
            gmeans_grp2.(nm1)=cmean;
            nm2=['std_' nm];
            gmeans_grp2.(nm2)=cstd;
            clear cmean nm nm1 nm2
        end
        clear d Md
    end
end

if summaryonly==0
%ANOVAN
factors=cell.empty;
fnames=cell.empty;
if grpF == 2
    factors=vertcat(factors,{allgrp});
    fnames=vertcat(fnames,{'group'});
end
if rwF > 1
    factors=vertcat(factors,{allrw});
    fnames=vertcat(fnames,{'row'});
end
if clF > 1
    factors=vertcat(factors,{allcl});
    fnames=vertcat(fnames,{'col'});
end

for cj=1:length(ciu)
    if length(factors) == 1
        [p(:,cj),tbl{1,cj},stats{1,cj}]=anovan(alldata{cj,1},factors,'varnames',fnames);
    elseif length(factors)>1
        [p(:,cj),tbl{1,cj},stats{1,cj}]=anovan(alldata{cj,1},factors,'model','interaction','varnames',fnames);
    end
    tblmodel=cell2table(tbl{1,cj});
    writetable(tblmodel, [outprefix '_module' num2str(cj) '_anovan'],FileType="spreadsheet")
    clear tblmodel
    close all hidden
    % multiple comparisons (only if main effects are present)
    dim = find(p(1:length(factors),cj)<0.05);
    if ~isempty(dim)
        [mcres{1,cj},~,f,gnames{1,cj}] = multcompare(stats{1,cj},'Dimension',dim);
        tblmc = array2table(mcres{1,cj},'VariableNames', ...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
        tblmc.("Group A")=gnames{1,cj}(tblmc.("Group A"));
        tblmc.("Group B")=gnames{1,cj}(tblmc.("Group B"));
        if exist('outprefix','var')
            fileout = [outprefix '_module' num2str(cj) '_multcomp_table'];
        else
            fileout = ['refMod_comsummary_' '_module' num2str(cj) 'multcomp_table'];
        end
        ver = length(dir(fullfile(pwd,[fileout '*'])));
        if ver == 0
            writetable(tblmc,fileout,'FileType','spreadsheet')
        else
            writetable(tblmc,[fileout '_run' num2str(ver+1),'.txt'],'FileType','spreadsheet')
        end
        clear ver
    end
end
end
%
% save the means and stdevs
if exist('outprefix','var')
    fileout = [outprefix '_comsummary'];
else
    fileout = 'refMod_comsummary';
end

ver = length(dir(fullfile(pwd,[fileout '*'])));
if ver == 0
    writetable(gmeans,fileout,'FileType','spreadsheet')
else
    writetable(gmeans,[fileout '_run' num2str(ver+1)],'FileType','spreadsheet')
end
clear ver
%% Write data for group2
if grpF == 2
    if exist('outprefix','var')
        fileout = [outprefix '_grp2_comsummary'];
    else
        fileout = 'refMod_grp2_comsummary';
    end

    ver = length(dir(fullfile(pwd,[fileout '*'])));
    if ver == 0
        writetable(gmeans_grp2,fileout,'FileType','spreadsheet')
    else
        writetable(gmeans_grp2,[fileout '_run' num2str(ver+1)],'FileType','spreadsheet')
    end
    clear ver
end










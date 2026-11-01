dirs = {'hemibrain',    'maleCNS_L',    'maleCNS_R'};

minOverlap = Inf;
maxOverlap = -Inf;
allOverlaps = zeros(length(dirs),3,3,2);
overlapsMBON01vsMBON12 = zeros(length(dirs),2);
for d = 1:length(dirs)
    % read the 3 MBON connectivity tables for this folder
    files = dir(fullfile(dirs{d},'*_KCinputs.csv'));
    numFiles = length(files);
    data = cell(numFiles,1);
    for i=1:numFiles
        data{i} = readtable(fullfile(files(i).folder,files(i).name)); %,optsA);
        data{i} = sortrows(data{i},'id');
    end

    % combine the two MBON12 files
    % this is kind of a hack, it would be better to do it by comparing the
    % names. It assumes the files are in alphabetical order and they are
    % the only 3 files in the folder matching *_KCinputs.csv
    MBON01 = data{1};
    MBON12 = [data{2}; data{3}];
    MBON12 = groupsummary(MBON12,{'id','type'},'sum','num_connections');
    MBON12.Properties.VariableNames{end} = 'num_connections';

    % find KC inputs that are in common between MBON01 and MBON12
    [c,ia,ib] = intersect(MBON01.id, MBON12.id);
    % figure('Name',dirs{d})
    % scatter(MBON01.num_connections(ia), MBON12.num_connections(ib), 'MarkerEdgeAlpha',0.1);
    overlaps = [numel(c)/height(MBON01), numel(c)/height(MBON12)];
    overlapsMBON01vsMBON12(d,:) = overlaps;

end

disp('rows: hemibrain, maleCNS_L, maleCNS_R')
disp('cols: as fraction of MBON01, as fraction of MBON12')
disp(overlapsMBON01vsMBON12)

%M = csvread("csv/0_101.csv")
%M2 = csvread("csv/9_69.csv")


%adjacency matrix AM containing all Graphs 
% 1 graph per cell
files = dir('*.csv');
L = length (files);
AM = cell(1,L);
filn = fopen('names_idx.txt', 'a+');

for i=1:L
    baseFileName = files(i).name;
    AM{i} = csvread(files(i).name);
    fprintf(filn, '%s\n', baseFileName)
    
end
% 
% files = dir('*.csv');
% for file = files'
%     csv = load(file.name)
%     %fprintf(filn, '%s\n', csv) 
%     % Do some stuff
% end

fclose(filn);
AM = transpose(AM);



new= cell2struct(AM(:,1)',{'am'});


save('CarboDimersVdwP2.mat', 'new', '-v7.3');
% filn2='AM_sparse3.mat';
% save filn2;
function [MS,S,attributes] = ExtractRiverNetwork(DEM,Amin)
%ExtractRiverNetwork extracts the river network of a sub-basin based on the
%sub_basin elevations (DEM) and computes attribute table of reaches (river
%segments between two consecutive confluences)
%
% Last update: 31/03/2016 (GC)
%
% Inputs:
% 1) DEM (GRIDobj)
% 2) Amin (number of cells) minimum drainage area
%
% Outputs:
% 1) MS (map struct) can be useful to obtain shapefilw
% 2) S (STREAMobj)
% 3) attributes -> matrix of reaches attributes:
%       columns of attributes matrix:
%        1 - Strahler order
%        2 - reach id
%        3 - id FN
%        4 - id TN
%        5 - x_FN
%        6 - y_FN
%        7 - elevation_FN      elevation(FN_indexes(:,1))
%        8 - x_TN
%        9 - y_TN
%       10 - elevation_TN      elevation(TN_indexes(:,1))
%       11 - drainage_area     (drainage_area(FN_indexes(:,1))+drainage_area(TN_indexes(:,1)))/2
%       12 - lengths
%       13 - slopes


%% Preprocessing

%DEM = fillsinks(DEM); %useless: included in option carve of function
%FLOWobj

%% Flow Directions

FD = FLOWobj(DEM,'preprocess','carve'); %includes flats preprocessing

%% Flow Accumulation

A  = flowacc(FD); %flow accumulation matrix
W = A>Amin;

%% Stream Network (reticolo idrografico)

S = STREAMobj(FD,W);
S = klargestconncomps(S,1); %keep only the main stream network

%% Compute Strahler order

 s = streamorder(S);
 % figure; [h,hs] = streamorder(S,'plot');

%% Obtain MS (31-3-16 WITH STRAHLER ORDER)

%impose reaches approximate length
reach_length = 10000000000000000000000; %very high length just to set one
MS = STREAMobj2mapstruct(S,'seglength',reach_length,...
'attributes',{'Strahler_order' s @median});
%                           'drainage_area' DA @mean});    %potrei
%                           calcolarla anche così (NB: prima però devo
%                           avere la matrice DA

%% ev. calcolare concavity index e ksn values 
%(vedi userguide 3 da riga 78)

%% Attribute table

% Syntax: [x,y,varargout] = STREAMobj2XY(S,varargin)
G = gradient8(DEM); %slope matrix

R = 6371; %mean earth radius (km)
cell_length = R*FD.cellsize*pi/180; %lunghezza del lato della cella (cellsize è in gradi sessagesimali)
DA = A.*(cell_length^2); %drainage area (km2)
%(A = number of upstream cells, to be multiplied per the dimension of each cell)

[x,y,elevation,drainage_area,slope] = STREAMobj2XY(S,DEM,DA,G);
%NB: la slope così calcolata si riferisce al singolo punto, non la uso

%dai vettori appena calcolati poi dovrei poter estrarre i valori in
%corrispondenza dei nodi, che sono nella tabella di S (in cui ci sono le
%coordinate x e y che mi permetteranno di riconoscerli)

%%

reach_ID = [1:length(MS)]';

%Find From Node (FN)

x_FN = nan(length(MS),1);
for i=1:length(MS)
    x_FN(i) = MS(i).X(1);
end

y_FN = nan(length(MS),1);
for i=1:length(MS)
    y_FN(i) = MS(i).Y(1);
end

%Find To Node (TN)

x_TN = nan(length(MS),1);
for i=1:length(MS)
    x_TN(i) = MS(i).X(end); %l'ultimo non è più nan -> prendo l'ultimo
end

y_TN = nan(length(MS),1);
for i=1:length(MS)
    y_TN(i) = MS(i).Y(end); %l'ultimo non è più nan -> prendo l'ultimo
end

%% Compute nodes indexes

%indexes of FN and TN nodes in the arrays x,y,elevation,drainage_area
%obtained through function STREAMobj2XY

FN_indexes = nan(length(x_FN),2);
for i=1:length(x_FN)
    indici = find( and( x==x_FN(i) , y==y_FN(i) ) );
    FN_indexes(i,1) = indici(1);
    FN_indexes(i,2) = indici(end); 
end

TN_indexes = nan(length(x_TN),2);
for i=1:length(x_TN)
    indici = find( and( x==x_TN(i) , y==y_TN(i) ) );
    TN_indexes(i,1) = indici(1);
    TN_indexes(i,2) = indici(end); 
end

%% Compute attributes matrix B
B = [reach_ID, x_FN, y_FN, elevation(FN_indexes(:,1)) , ...
    x_TN, y_TN, elevation(TN_indexes(:,1)) ,...
    (drainage_area(FN_indexes(:,1))+drainage_area(FN_indexes(:,1)))/2]; 
%NB: nomore assume that drainage area of a 
%reach is the drainage area of its upstream node (FN) -> 31-3-16: take the
%mean of FN and TN values

%mancano lunghezze e pendenze

%ricavo gli indici dei nodi in S confrontando le coordinate

FN_indexes_S = nan(length(x_FN),2);
for i=1:length(x_FN)
    indici = find( and( S.x==x_FN(i) , S.y==y_FN(i) ) );
    FN_indexes_S(i,1) = indici(1);
    FN_indexes_S(i,2) = indici(end); %invece che lasciare nan riscrive l'indice uguale
end

TN_indexes_S = nan(length(x_TN),2);
for i=1:length(x_TN)
    indici = find( and( S.x==x_TN(i) , S.y==y_TN(i) ) );
    TN_indexes_S(i,1) = indici(1);
    TN_indexes_S(i,2) = indici(end); %invece che lasciare nan riscrive l'indice uguale
end

lengths = S.distance(FN_indexes_S(:,1)) - S.distance(TN_indexes_S(:,1));

slopes = (elevation(FN_indexes(:,1)) - elevation(TN_indexes(:,1)))./lengths;

B = [B lengths slopes];
% columns of B matrix:
% 1-reach_ID
% 2-x_FN
% 3-y_FN
% 4-elevation_FN      elevation(FN_indexes(:,1))
% 5-x_TN
% 6-y_TN
% 7-elevation_TN      elevation(TN_indexes(:,1))
% 8-drainage_area     
% 9-lengths
% 10-slopes

%% check directions

%wrong_dir = find( elevation(FN_indexes(:,1)) < elevation(TN_indexes(:,2)) );

%% FN-TN matrix

uniqueNodes = unique( [ x_FN,y_FN ; x_TN,y_TN ] , 'rows'); 

idNodes = [1:length(uniqueNodes)]';

uniqueNodes = [idNodes, uniqueNodes];

idFN = nan(length(x_FN),1);
for i=1:length(x_FN)
    FN_index = find(  and( uniqueNodes(:,2)==x_FN(i) , uniqueNodes(:,3)==y_FN(i) )  );    
    idFN(i)= uniqueNodes(FN_index,1);
    %sarebbe stato equivalente fare così
    %idFN(i)= find(  and( uniqueNodes(:,2)==x_FN(i) , uniqueNodes(:,3)==y_FN(i) )  );
    %visto che idNodes coincide con l'indice
end

idTN = nan(length(x_TN),1);
for i=1:length(x_TN)
    TN_index = find(  and( uniqueNodes(:,2)==x_TN(i) , uniqueNodes(:,3)==y_TN(i) )  );
    idTN(i)= uniqueNodes(TN_index,1);
end

FN_TN_messy_matrix = [ reach_ID idFN idTN ];


%% reassign node IDs (riordino nodi)

[ new_idFN, new_idTN, outlet_node_id ] = reassignNodeIDs(idFN, idTN );

FN_TN_ordered_matrix = [ reach_ID new_idFN new_idTN ]; 

new_uniqueNodes = unique( [ new_idFN, x_FN, y_FN ; ...
    new_idTN, x_TN, y_TN ] , 'rows'); %give new ids to coordinates couples


%% new attribute matrix
%NB non dovrebbe essere cambiato l'ordine degli archi, ma solo il nome dei
%nodi

attributes = [FN_TN_ordered_matrix B(:,2:end)];

% columns of attributes matrix:
% 1-reach id
% 2-id FN
% 3-id TN
% 4-x_FN
% 5-y_FN
% 6-elevation_FN      elevation(FN_indexes(:,1))
% 7-x_TN
% 8-y_TN
% 9-elevation_TN      elevation(TN_indexes(:,1))
% 10-drainage_area     
% 11-lengths
% 12-slopes

%% add attributes to MS

for i=1:length(MS)
    MS(i).reach_id = attributes(i,1);
    MS(i).id_FN = attributes(i,2);
    MS(i).id_TN = attributes(i,3);
    MS(i).x_FN = attributes(i,4);
    MS(i).y_FN = attributes(i,5);
    MS(i).elevationFN = attributes(i,6);
    MS(i).x_TN = attributes(i,7);
    MS(i).y_TN = attributes(i,8);
    MS(i).elevationTN = attributes(i,9);
    MS(i).drainage_area = attributes(i,10);
    MS(i).length = attributes(i,11);
    MS(i).slope = attributes(i,12);
end

%% add strahler order to attributes
N = length(MS);
strahler = nan(N,1);
for i = 1:N
    strahler(i) = MS(i).Strahler_order;
end

attributes = [strahler attributes];

% columns of new attributes matrix:
%        1 - Strahler order
%        2 - reach id
%        3 - id FN
%        4 - id TN
%        5 - x_FN
%        6 - y_FN
%        7 - elevation_FN      elevation(FN_indexes(:,1))
%        8 - x_TN
%        9 - y_TN
%       10 - elevation_TN      elevation(TN_indexes(:,1))
%       11 - drainage_area     (drainage_area(FN_indexes(:,1))+drainage_area(TN_indexes(:,1)))/2
%       12 - lengths
%       13 - slopes


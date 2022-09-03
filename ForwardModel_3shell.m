function [gridpoints,TransferMat] = ForwardModel_3shell(Resolution, ModelParams)

%Grid Locations
Radius = ModelParams.R(1) ;
[X,Y,Z] = meshgrid(-Radius:Resolution:Radius,-Radius:Resolution:Radius,-Radius*.1:Resolution:Radius) ;
X = reshape(X,1,[]) ;
Y = reshape(Y,1,[]) ;
Z = reshape(Z,1,[]) ;
grid = find(X.^2+Y.^2+Z.^2<Radius^2) ;
gridpoints = [X(grid); Y(grid); Z(grid)] ;
clear('X','Y','Z') ;
GridLocs_temp = reshape(gridpoints,1,3,[]) ;
GridLocs_temp = repmat(GridLocs_temp,[21,1,1]) ;

%Electrode Positions
load('ElecPosXYZ') ;
ElectrodePos = [] ;
for i=1:21
    ElectrodePos(i,:) = Radius*ElecPos{i}.XYZ ;
end
ElectrodePos = repmat(ElectrodePos,[1,1,size(GridLocs_temp,3)]) ;

TransferMat = zeros(size(GridLocs_temp)) ;

for shell=1:3
    GridLocs = ModelParams.Mu(shell)*GridLocs_temp ;
    
    C1 = 2*sum((ElectrodePos-GridLocs).*GridLocs,2)./(sum((ElectrodePos-GridLocs).^2,2)).^1.5 + ...
    1./sqrt(sum((ElectrodePos-GridLocs).^2,2))- ...
    1./sqrt(sum((ElectrodePos).^2,2)) ;
    C2 = 2./(sum((ElectrodePos-GridLocs).^2,2)).^1.5 + ...
    (sqrt(sum((ElectrodePos-GridLocs).^2,2))+ sqrt(sum((ElectrodePos).^2,2)))./...
    (sqrt(sum((ElectrodePos).^2,2)).*sqrt(sum((ElectrodePos-GridLocs).^2,2)).*...
    (sqrt(sum((ElectrodePos).^2,2)).*sqrt(sum((ElectrodePos-GridLocs).^2,2))+...
    sum((ElectrodePos).^2,2)-sum(GridLocs.*ElectrodePos,2))) ;
    h = (repmat(C1-C2.*sum(GridLocs.*ElectrodePos,2),[1,3,1]).*GridLocs + ...
        repmat(C2.*sum(GridLocs.^2,2),[1,3,1]).*ElectrodePos)./...
        (repmat(4*pi*ModelParams.Sigma(3)*sum(GridLocs.^2,2),[1,3,1])) ;
    TransferMat = TransferMat + ModelParams.Lambda(shell)*h ;
end

TransferMat = reshape(TransferMat,21,[]) ;

    




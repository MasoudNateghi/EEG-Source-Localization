function Display_Potential_3D(Radius,ElecPot)

load('ElecPosXYZ') ;
load('ElecPatch') ;

%Nasion, Inion, Ears
Nz = Radius*[0.927+.25	 -0	-0.375] ;
Iz = Radius*[-0.906-.4	-1.11e-16 -0.423] ;

Ver = [] ;
for i=1:21
    Ver(i,:) = Radius*ElecPos{i}.XYZ ;
end
Face = ElecPatch ;

Cmat = reshape(ElecPot,[],1) ;
patch('Vertices',Ver,'Faces',Face,'FaceVertexCData',Cmat,'FaceColor','interp','FaceAlpha',1)
hold on
for i=1:21
    text(Ver(i,1) ,Ver(i,2), Ver(i,3), ElecPos{i}.Name) ;
end
text(Nz(1),Nz(2),Nz(3), 'Nasion','HorizontalAlignment','center','EdgeColor','r') ;
text(Iz(1),Iz(2),Iz(3), 'Inion','HorizontalAlignment','center','EdgeColor','r') ;
xlim([-1.5,1.5]*Radius) ;
ylim([-1.2,1.2]*Radius) ;
zlim([-.5,1.3]*Radius) ;
view(90,90) ;
colorbar



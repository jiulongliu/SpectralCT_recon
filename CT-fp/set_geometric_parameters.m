function para=set_geometric_parameters(nx,ny,nt,nv,nd,offset,version,GPU)

if nargin < 8
    GPU = 1;
    if nargin < 7
        version = 1;
        if nargin < 6
            offset=0;
            if nargin<5
                nd = 512;
                if nargin<4
                    nv=360;
                    if nargin<3
                        nt=1;
                        if nargin<2
                            ny=256;
                            if nargin<1
                                nx=256;
                            end
                        end
                    end
                end
            end
        end
    end
end



% fan-beam geometry %
% nx=256;ny=nx;% number of pixels: nx*ny
dx=single(250/nx);dy=dx;% pixel size: dx*dy
% nv=668;% number of views (projections)

SO=single(1000);% the distance from source to isocenter
OD=single(500);% the distance from isocenter to detector
% nd=512;% number of detector
dy_det=single(0.388*1024/nd);% size of detector
sd_phi=single(2*pi/nv*(0:nv-1));% view angles
y_os=single(offset*dy_det);% isocenter offset with respect to the detector center
y_det=single(((-nd/2:nd/2-1)+0.5)*dy_det)+y_os;% detector coordinate
y_det2=single((-nd/2:nd/2)*dy_det)+y_os;
% normalize so that dx=dy=1 %
scale=dx;
SO=SO/scale;OD=OD/scale;y_det=y_det/scale;
dy_det=dy_det/scale;y_os=y_os/scale;
% parameters (single image frame) %
% note: this version works for multiple images as well by supplying id_v -- which projection views are for each frame%
% load phantom_2d.mat x0 % 2D phantom image
% nt=1;
% X0=x0(:);
Id_v=cell(1,nt);
Id_v{1}=uint32(0:nv-1);

id_X=[];Nv=zeros(1,nt);
for i=1:nt
    id_X=[id_X (i-1)*ones(1,numel(Id_v{i}))];
    Nv(i)=numel(Id_v{i});
end
tmp_size=max(Nv);
id_Y=zeros(tmp_size,nt);
for i=1:nt
    id_Y(1:Nv(i),i)=Id_v{i};
end
para=struct('version',uint32(version),'GPU',uint32(GPU),...
    'SO',single(SO),'OD',single(OD),'scale',single(scale),'nx',uint32(nx),'ny',uint32(ny),'nv',uint32(nv),'nd',uint32(nd),...
    'sd_phi',sd_phi,'y_det',y_det,'id_X',uint32(id_X),'nt',uint32(nt),'dy_det',dy_det,'y_os',y_os,...
    'id_Y',uint32(id_Y),'Nv',uint32(Nv),'tmp_size',uint32(tmp_size),...
    'cos_phi',cos(sd_phi),'sin_phi',sin(sd_phi),'cos_det',[],'sin_det',[]);        
angle_det=atan2(y_det,SO+OD);para.cos_det=cos(angle_det);para.sin_det=sin(angle_det);
angle_det2=atan2(y_det2,SO+OD);para.cos_det2=single(cos(angle_det2));para.sin_det2=single(sin(angle_det2));% for finite-size beam


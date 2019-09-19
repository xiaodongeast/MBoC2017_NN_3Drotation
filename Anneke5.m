% need to modify the rotate and the sort so the output is in a customized
% field. And then i can use the same function for the random point
% generated
% 1 is cell; 2 is frag ;3 r, 4-6 xyz for frag 7 MT no, 8,9,10 MT x,y,z 





%
function [subT,rots2,ac] = Anneke5( D,select )
%% make the modification to disable all the useless ones, 
% leave 3 functions: random point generation; neareast neigbor calcualtion;
% selection the size of the golgi donuts.
%%
%for anneke's Golgi MT hotspot
%Based on size and position of the Golgi mini-stack, rescale and shift the
%coordinate.
% Project the position of MTs to the surface of the Golgi 
% For MTs on one Golgi fragment, find the MTs has the mini nearest neighbor distance[x,y,z]
% Calculate the rotation matrix between[0,0,1] and the point [x,y,z]; and use these to calculate the new position 
% [x_new; y_new; z_new]=[rotation matrix]*[x,y,z]

%    in the result structure rots2: 'mtRotate' is the final cooridnate of the
%    mt after rescale, shifting,rotation, 
%   'randRotate is the final coordinate of the random pint after rescale,
%   shifting and rotation
% the last column of the mtSort is the neareast distance for each point of
% the mt after shifting rescaling. point in the first row has the min of the nearest
% distance.
% similary 'randSort' is the sorting for the nearest neigbour

% %first use ogSplit to organize data into a structure based on the fragment
% No. which is the identification key

%%  
s=ogSplit(D);
% mtCor is the orignal data

%% based on the fragment r size, the coordinate is rescale and translate to
%the center (0,0,0);
s2=transform(s);
%nmtCor is the data after shifting the rescale
%% project to the surface 
s3=transProj(s2);

%nmtcor_proj is the new coordinate after projection
%%  calculate the mtNearest distance. Here the original corordinte was used, change to the projected one if need
ss=sortNear(s3,'mtCor','mtSort','mtNearest','mtIPDM');

%% generate random points, change the code in the function already. So this program is now independent of the modification above_03042017
rs=randPoint(ss);

%% from here all the program disabled to avoid the plot
%ss=1;


%%
%rots=rotation(ss,'mtSort','mtRotate');
%rots2=rotation2(rots,'mtRotate','mtRotate2');
% in the structure 'rots', I have finished all the change for the original
% data
%% calculate the distance for the random one
ss2=sortNear(rs,'randCtr','randSort','randNearest','randIPDM');

%% disable the rotation
%rots2=rotation(ss2,'randSort','randRotate');
%rots2=rotation2(rots2,'randRotate','randRotate2');
%% for the selection. I have to  add rots2=ss2, this in order to make the program run without other modification
rots2=ss2;    % to avoid rotation

%close all;
% with the condition
subT=rots2;
if select(1)>0
    subT=struct2table(rots2);
  subT=subT(subT.r>select(1),:);
  subT=subT(subT.r<select(2),:);
  subT=table2struct(subT);
end

%plotPoint(rots2,i); % use the whole to plot the select i;
%totalPoint=finalPlot(subT,'mtRotate2','randRotate2'); % use the slected for the further cal   
% Change to allow only 1 fragment be selected.
%[IPDM, Nearest]=calDis(subT);
%[h,p]=ttest(Nearest(:,1),Nearest(:,2));
ac=copyDs(subT,'mtNearest');
end
%% have a function organize the data,
% output as a structure, S.cell; S.frag;S.center,S.Index, S.MTposition

function [Ds]=ogSplit(D)
mtStartCol=8;
mtStopCol=10;
rStartCol=4;
rStopCol=6;
constantRow=1;
 Dt=array2table(D,'VariableNames',{'celln','fNo','r','rx','ry','rz','MTNo','mx','my','mz'});
 Ds=struct('celln',[],'fNo',[],'r',[],'rCor',[],'mtCor',[]);
[G,ID]=findgroups(Dt.fNo);
for i=1:length(ID)
    A=Dt(Dt.fNo==ID(i),:);
     a=table2array(A);
     am=a(:,mtStartCol:mtStopCol);
     ar=a(constantRow,rStartCol:rStopCol);
     Ds(i).rCor=ar;
     Ds(i).mtCor=am;
     Ds(i).r=A.r(constantRow);
     Ds(i).celln=A.celln(constantRow);
     Ds(i).fNo=A.fNo(constantRow);          
end
end


%% merge the stuctures back into a table 
function [Dn2]=ogMerge(D)
end
%%
%% transform1 add the new column gives resclale with the r, shift to the center, move the points to the surface 
%% based on the center rescale, and move the coordinate.
function [ds2]=transform(ds)

ds2=ds;
for i=1:length(ds)
    scale=ds(i).r;  % scale is based on distance r;
    
   ds2(i).nr=ds(i).r/scale;
    ds2(i).nrCor=ds2(i).rCor./scale
    ds2(i).nmtCor=ds2(i).mtCor./scale
    % scae all the cooridate, rcor, mt cor;
    % next shift;
    
      shift=ds2(i).nrCor;
    shift=shift.*-1;
    [m,n]=size(ds2(i).mtCor);
    shiftM=ones(m,1)*shift;
   ds2(i).shiftM=shiftM;
   ds2(i).nrCor=ds2(i).nrCor+shift;
    ds2(i).nmtCor=ds2(i).nmtCor+shiftM;
end   

end
%% project only the mtCoordinate to move it to the suface
function[ds2]=transProj(ds)
ds2=ds;
for i=1:length(ds)
    m= ds2(i).nmtCor;
    m=m';
    m2 = bsxfun(@rdivide, m, sqrt(sum(m.^2,1)));
    ds2(i).nmtCor_Proj=m2';
end
end
 

%% add the rotation to move the coordinate

function [ds2o]=rotation(dsi,colS,colR)
ds2o=dsi;
rotateZ=[0,0,1];
for iii=1:length(dsi)
   mtemp=dsi(iii).(colS);
   mtemp=mtemp(:,1:3)';
    r = vrrotvec(mtemp(:,1),rotateZ);
    rotateM=vrrotvec2mat(r);
    mtemp2=rotateM*mtemp;
    ds2o(iii).(colR)=mtemp2';
    % abvoe finish the rotation along the first point
   
end
end
%% add the function for the 2nd rotation. here Irina want to rotate again 
% and use the 2nd point as the refer. so the 2nd point is in one logitude.


function [dsr2]=rotation2(dsi,cols,colR)
dsr2=dsi;
for iij=1:length(dsi)
   a=dsi(iij).(cols);
   a=a';
 oldP=a(:,2);
 % compute the r ofthe circle at the given z.
     r=sqrt(1-oldP(3)^2);
     % the vector to rotate, I will not change z, so z is 0; I will change
     % the y value to the 0, and the x will remain unchanged. so restrict
     % the direction of the rotation
    refer=[r,0,0];
    % compute the rotation need
     ro = vrrotvec([oldP(1),oldP(2),0],refer);
     % covert to the regulor rotation matrix
    rotateM=vrrotvec2mat(ro);
     newP=rotateM*a;
     dsr2(iij).(colR)=newP';
end
end

%% 3D plot with a sphere

function []=plotPoint(test,i)

figure(1);
[x,y,z]=ellipsoid(0,0,0,1,1,1,20);
s=surf(x,y,z);
alpha(s,0.5);
s.FaceColor=[0.5,0.5,0.5];
hold;
scatter3(test(i).mtCor(:,1),test(i).mtCor(:,2),test(i).mtCor(:,3),'ro');
scatter3(test(i).mtCor(:,1),test(i).mtCor(:,2),test(i).mtCor(:,3),'b*');
[x1,y1,z1]=ellipsoid(test(i).rCor(1),test(i).rCor(2),test(i).rCor(3),test(i).r,test(i).r,test(i).r,20);
ss=surf(x1,y1,z1);
alpha(ss,0.5);
scatter3(test(i).nmtCor(:,1),test(i).nmtCor(:,2),test(i).nmtCor(:,3),'b*');
scatter3(test(i).nmtCor_Proj(:,1),test(i).nmtCor_Proj(:,2),test(i).nmtCor_Proj(:,3),'r*');
scatter3(test(i).randCtr(:,1),test(i).randCtr(:,2),test(i).randCtr(:,3),'g*');
 scatter3(test(i).mtSort(1,1),test(i).mtSort(1,2),test(i).mtSort(1,3),100,'d','MarkerFaceColor',[1 1 .75],'MarkerEdgeColor','k');
figure(2);
[x,y,z]=ellipsoid(0,0,0,1,1,1,20);
s=surf(x,y,z);
alpha(s,0.5);
s.FaceColor=[0.5,0.5,0.5];
hold;
scatter3(test(i).mtRotate(:,1),test(i).mtRotate(:,2),test(i).mtRotate(:,3),100,'d','MarkerFaceColor',[1 1 .75],'MarkerEdgeColor','k')
scatter3(test(i).randRotate(:,1),test(i).randRotate(:,2),test(i).randRotate(:,3),100,'*','MarkerFaceColor',[1,0,1],'MarkerEdgeColor',[1,0,0])
figure(3);
[x,y,z]=ellipsoid(0,0,0,1,1,1,20);
s=surf(x,y,z);
alpha(s,0.5);
s.FaceColor=[0.5,0.5,0.5];
hold;
scatter3(test(i).mtRotate2(:,1),test(i).mtRotate2(:,2),test(i).mtRotate2(:,3),100,'d','MarkerFaceColor',[1 1 .75],'MarkerEdgeColor','k')
end
%% random point generate
%gnerate random point
% and fill to a new structure
function [totalPoint]=finalPlot(ds,col,col2)
% change this to allow only one fragment was selected.
figN=4;
figure(figN);
hold;
title('mt relative position--- from experiment');
[x,y,z]=ellipsoid(0,0,0,1,1,1,20);
s=surf(x,y,z);
alpha(s,0.5);
%s.FaceColor=[0.5,0.5,0.5];
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis, imageHeight');
if length(ds)>1
sss=struct2table(ds);
randr=sss.(col2);
mtr=sss.(col);
mt=cell2mat(mtr);
randr=cell2mat(randr);
else sss=ds;
    mt=sss.(col);
    randr=sss.(col2);
end
%randr=sss.(col2);
% this should be only when more than 2 were selected.
% mtr,mt are values need to be handled when there is only one was selected
%mtr=sss.(col);
 

%mt=cell2mat(mtr);

[mtPointN,whatever]=size(mt);
scatter3(mt(:,1),mt(:,2),mt(:,3),'d','MarkerFaceColor','none','MarkerEdgeColor',[1,0,1],'LineWidth',1.5)

figure(figN+1);
hold;
title('random position---mataching microtubule relative position from experiment');
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis, imageHeight');
[x,y,z]=ellipsoid(0,0,0,1,1,1,20);
s=surf(x,y,z);
alpha(s,0.5);
%s.FaceColor=[0.5,0.5,0.5];

%randr=sss.(col2);
%randr=cell2mat(randr);
[randPointN,whatever]=size(randr);
scatter3(randr(:,1),randr(:,2),randr(:,3),'d','MarkerFaceColor','none','MarkerEdgeColor',[0.8,0.2,0],'LineWidth',1.5)
[totalPoint]=[mtPointN,randPointN];
end

function[ipdm, near]= calDis(ds)
if length(ds)>1    
sss=struct2table(ds);
mtIPDM=sss.mtIPDM;
mtIPDM=cell2mat(mtIPDM);
randIPDM=sss.randIPDM;
randIPDM=cell2mat(randIPDM);
mtNear=sss.mtNearest;
mtNear=cell2mat(mtNear);
randNear=sss.randNearest;
randNear=cell2mat(randNear);
else
    sss=ds;
    mtIPDM=sss.mtIPDM;
      randIPDM=sss.randIPDM;
      mtNear=sss.mtNearest;
      randNear=sss.randNearest;
      
    
end


ipdm=[mtIPDM,randIPDM];
near=[mtNear, randNear];
end


function[dsr]=randPoint(ds)
dsr=ds;
for i=1:length(ds);
    mt=ds(i).mtCor; % based on the number of the original size of MT number
    [m,n]=size(mt');
    r1=randn(m,n);
    r2 = bsxfun(@rdivide, r1, sqrt(sum(r1.^2,1)));
    r2=r2';
    dsr(i).randCtr=r2;
end

end

%%Sort the data according to the nearest neigbor distance
 

%%Sort the data according to the nearest neigbor distance

function[ds3,dDmin,dDu]=sortNear(rs1,colS,colR,colN,colIPDM)
    ds3=rs1;
    
  for ii=1:length(ds3);
    mmmm=ds3(ii).(colS);
    Index=4;
    dD=pdist2(mmmm,mmmm);
    IPDMu=triu(dD);
   dDu=reshape(IPDMu',numel(IPDMu),1);
   dDu=nonzeros(dDu);
  [im,in]=size(mmmm);
  dDmin=ones(im,1);
      for j=1:im;
             dDmin(j)=min(nonzeros(dD(j,:)));
      end
    
      dD2=[mmmm,dDmin];
     dD2=sortrows(dD2,Index);
     ds3(ii).(colR)=dD2;
     ds3(ii).(colN)=dDmin;
     ds3(ii).(colIPDM)=dDu;
  end
   
end



function[aa]= copyDs(ds,colR)
if length(ds)>1    
sss=struct2table(ds);
aa=sss.(colR);
aa=cell2mat(aa);
 
else
    sss=ds;
    aa=sss.mtIPDM;
      
    
end


 
end
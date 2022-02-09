function mothers=phy_setObjectLinks(object,channel,noiteration)
global segmentation

% new pedigree construction based on budnecks detection without using
% budneck mapping


% pedigree start frame can be different from movie start frame
%

%displayImage=segmentation.realImage;

phy_progressbar;

% firstMCells

firstFrame=segmentation.pedigree.start;

cells=segmentation.(object)(firstFrame,:);
pix=find([cells.ox]~=0);
n=[cells.n];


firstMCell=n(pix);
firstCells=[];

minDivisionTime=20;%segmentation.pedigree.minDivisionTime;

exclude=firstMCell;

tcells=segmentation.(['t' object])  ;
 
% init parentage -
for i=1:numel(tcells)
    if tcells(i).N~=0
        tcells(i).setMother(0);
    end
end




% sort tcells according to appearance timing;

order=[];
for i=1:numel(tcells)
    if numel(find(firstMCell==tcells(i).N))
           continue
    end
    
    if tcells(i).N~=0
        if tcells(i).mother==0
            if tcells(i).detectionFrame>=segmentation.pedigree.start && tcells(i).detectionFrame<segmentation.pedigree.end
                
               if ~numel(find(find(segmentation.discardImage)==tcells(i).detectionFrame))
                 
                pix=find(exclude==tcells(i).N);
                if  numel(pix)==0
                    order=[order; [i tcells(i).detectionFrame tcells(i).N]];
                end
               end
            end
        end
    end
    % i,a=order(i)
end

[narr sortindex]=sortrows(order,2);

phy_progressbar;
pause(0.1);

% first identify all possible candidates and rank them

%if numel(candarrstore)==0
    
    candarr=zeros(length(narr(:,1)),10);
    scorearr=zeros(length(narr(:,1)),10);
    
    for k=1:numel(narr(:,1))
        phy_progressbar(double(k)/numel(narr(:,1)));
        
        cindex=narr(k,1);
        
        fr=tcells(cindex).detectionFrame;
        
        cells1=segmentation.(object)(fr,:);
        ncells1=[cells1.n];
        pixcells1=find(ncells1>0);
        cells1=cells1(pixcells1);
        
        targetCell=tcells(cindex).Obj(1);
        score=zeros(1,max(narr(:,1)));
        
        [candidates dist]=findNeighbors(targetCell,cells1,50);
        
        if numel(candidates)==0
            continue
        end
        
        score(candidates)=1./dist';
        
        
        if channel~=0
        candidates=scoreDiv(tcells,candidates,fr,channel,minDivisionTime);
        
        score(candidates)=score(candidates)+1;
        end
        
        [ord ind]=sort(score,'descend');
       % ord
        
        pix=find(ord==0,1,'first');
        
        candarr(k,1:pix-1) =ind(1:pix-1);
        scorearr(k,1:pix-1)=ord(1:pix-1);
        
    end
    
    phy_progressbar(1);
    pause(0.1);

problems=1;
energy=0;
cc=1;

%[narr(:,3) candarr(:,1:3) scorearr(:,1:3)]

%return;

% detect and fix problems based on timings

listDau=[];

firstFrame=find(segmentation.([object 'Mapped']),1,'first'); %remove -1 ?
firstframe=segmentation.pedigree.start;

tic;

ener=[];
moth.mothers=[];
cc2=0;

ntcells=[tcells.N];

if nargin==3 % no iterations wanted : just build the pedigree tree based on scoring
    
    mothers=buildTree(narr,candarr,tcells,ntcells);

    [problems energy]=computeEnergy(mothers,minDivisionTime,10000,firstframe);
    

else
    
    cc=1;
    mothers=buildTree(narr,candarr,tcells,ntcells);
    [problems energy]=computeEnergy(mothers,minDivisionTime,10000,firstframe);
   % pause
    ener(cc)=energy;
    moth(cc).mothers=mothers;
    cc=2;
    
    fprintf('Now performing optimization....\n');
    
    while cc<10*length(narr(:,3)) 

        move=1;
        cc3=0;
        while move
        ind=randi(numel(problems));
        
        pix=find(narr(:,3)==problems(ind));
            
        ncandidates=find(candarr(pix,:));
        if length(ncandidates)>1
            break;
        end
        
        if cc3==1000
           break 
        end
        cc3=cc3+1;
        end
        
        if cc3==1000
            break;
        end
        
        
        listcand=candarr(pix,ncandidates);
        listcand=circshift(listcand',-1);
        listcand=listcand';
        
     %   narr(pix,3)
        
        %a=candarr(pix,:)
        
        candarr(pix,1:length(listcand))=listcand;
        
        %a=candarr(pix,:)
        %pause
        
        mothers=buildTree(narr,candarr,tcells,ntcells);
           
        [problems energy]=computeEnergy(mothers,minDivisionTime,10000,firstframe);
        %problems
        
       
            if energy<ener(cc-1)
            ener(cc)=energy;
            moth(cc).mothers=mothers;
            else
            q=exp(-(energy-ener(cc-1)/0.1));
            
            if rand(1)<q
                ener(cc)=energy;
                moth(cc).mothers=mothers;
            else
                ener(cc)=ener(cc-1);
                moth(cc).mothers=moth(cc-1).mothers;
            end
            
            end
            
           % a=ener(cc)
   
        cc=cc+1;
    end
    toc;
    
    
    
   
        %figure, plot(ener);
        [ccm cci]=min(ener);
        idx=cci;
        mothers=moth(idx).mothers;
        
        %problems
end

        
for i=1:numel(mothers)
    n=[tcells.N];
    list=mothers(i).daughterList;
    detect=mothers(i).budTimes;
    %tcells(i).setMother(0);
    
    for j=1:numel(list)
        dau=list(j);
        %n(i);
        
        if numel(find(tcells(i).daughterList==dau))==0
            pix=find(n==dau);
            tcells(pix).setMother(n(i));
            tcells(i).addDaughter(dau,tcells(pix).detectionFrame,tcells(pix).detectionFrame); %add a new daughter to the mother
        end
    end
end

fprintf('Done !\n');

%b=tcells(1)

function [candarr scorearr]=makeNewConfig(narrin,candarrin)

candarr=candarrin;
scorearr=scorearrin;

%problems,narrin(:,1)

dau2=find(narrin(:,3)==problems(2));
dau3=find(narrin(:,3)==problems(3));

temp2=candarr(dau2,:);
tempscore2=scorearr(dau2,:);

pix2=find(temp2~=0);
temp2=temp2(pix2);
tempscore2=tempscore2(pix2);

temp2=circshift(temp2',-1);
tempscore2=circshift(tempscore2',-1);

temp2=temp2';
tempscore2=tempscore2';


temp3=candarr(dau3,:);
tempscore3=scorearr(dau3,:);

pix3=find(temp3~=0);
temp3=temp3(pix3);
tempscore3=tempscore3(pix3);

temp3=circshift(temp3',-1);
tempscore3=circshift(tempscore3',-1);

temp3=temp3';
tempscore3=tempscore3';


r= tempscore3(1)/(tempscore3(1)+tempscore2(1));

if rand(1)<r
    candarr(dau3,pix3)=temp3;
    scorearr(dau3,pix3)=tempscore3;
else
    candarr(dau2,pix2)=temp2;
    scorearr(dau2,pix2)=tempscore2;
end


function  mothers=buildTree(narr,candarr,tcells,ntcells)

mothers.daughterList=[];
mothers.budTimes=[];
mothers.detectionFrame=[];
mothers.n=[];

for i=1:numel(tcells)
    mothers(i).detectionFrame=tcells(i).detectionFrame;
    mothers(i).n=tcells(i).N;
end

n=ntcells;

for i=1:length(narr(:,1))
    if candarr(i,1)~=0
        ind=find(n==candarr(i,1));
        
        mothers(ind).daughterList=[mothers(ind).daughterList narr(i,3)];
        mothers(ind).budTimes=[mothers(ind).budTimes narr(i,2)];
    end
end

function [out distout]=findNeighbors(targetCell,cellsin,dist)

mx=[cellsin.ox];
mx=repmat(mx',[1 size(mx,2)]);
mx=mx-mx';

my=[cellsin.oy];
my=repmat(my',[1 size(my,2)]);
my=my-my';

sz=sqrt(mean([cellsin.area]));

d=sqrt(mx.^2+my.^2);
pix=d<dist;
pix=pix & tril(ones(size(d)),-1);

[row,col] = find(pix);

if numel(pix)==0
    out=[];
    return
end

nc=[cellsin.n];


val= find(nc==targetCell.n);


pix=[find(col==val) ; find(row==val)];

col=col(pix);
row=row(pix);

n=length(row);
%
fuse=[];
dist=[];
%


% find min distances between cells
for i=1:n
    
    
    x1=cellsin(row(i)).x;
    if size(x1,1)~=1
    x1=x1';
    end
        
    x2=cellsin(col(i)).x;
    if size(x2,1)~=1
    x2=x2';
    end
    
    y1=cellsin(row(i)).y;
    if size(y1,1)~=1
    y1=y1';
    end
    
    y2=cellsin(col(i)).y;
    if size(y2,1)~=1
    y2=y2';
    end
    
    %row(i)
    % line(cellsin(row(i)).x,cellsin(row(i)).y,'Color','r','Marker','o');
    
    x1p=repmat(x1',[1 size(x1,2)]);
    x2p=repmat(x2',[1 size(x2,2)]);
     %i,size(x1p),size(x2p)

%     if size(x1p)==0 | size(x2p)==0
%        out=[];
%        distout=[];
%        return
%     end
        x=x1p-x2p';
    y1p=repmat(y1',[1 size(y1,2)]);
    y2p=repmat(y2',[1 size(y2,2)]);
    y=y1p-y2p';
    
    d=sqrt(x.^2+y.^2);
    %
    %row(i),min(min(d))
    
    %pix=d<20;
    %pix=pix & ~diag(ones(1,size(d,1))) ;%tril(ones(size(d)),-1);
    
    %pix=find(pix);
    
   % if numel(pix)>0
        %fuse=[fuse i];
        dist=[dist min(d(:))];
   % end
end


out=[];
distout=[];


    for i=1:n
    %val,row(fuse(1))
    if row(i)==val
        out=[out ; col(i)];
        distout=[distout; dist(i)];
    else
        out=[out ; row(i)];
        distout=[distout; dist(i)];
    end
    end
    

 %   out
%out=unique(out)

    out=nc(out);

function out=scoreAxis(cells1,candidates,target)

out=[];

xarr=target.x;
yarr=target.y;

sizex=round(max(xarr)-min(xarr)+50);
sizey=round(max(yarr)-min(yarr)+50);

BW=poly2mask(xarr-min(xarr)+25,yarr-min(yarr)+25,sizex,sizey);

stat=regionprops(BW,'Eccentricity','Orientation');

if stat.Eccentricity>0.25
    dist=30;
    
    vec=-dist:5:dist;
    linex=mean(xarr)+vec*cos(-stat.Orientation*2*pi/360);
    liney=mean(yarr)+vec*sin(-stat.Orientation*2*pi/360);
    
    %figure, imshow(BW); hold on;
    %line(linex-min(xarr)+25,liney-min(yarr)+25,'Color','m');
    
    %figure;
    
    for i=1:length(cells1)
        % candidates,cells1(i).n
        pix=find(candidates==cells1(i).n);
        if numel(pix)~=0
            % pix
            %   line(cells1(i).x,cells1(i).y); hold on;
            
            
            if mean(inpolygon(linex,liney,cells1(i).x,cells1(i).y))~=0
                out=[out cells1(i).n];
            end
        end
    end
    
    %line(xarr,yarr,'Color','g'); hold on
    %        line(linex,liney,'Color','r'); hold on
    %title(num2str(target.n))
    %axis square
    
end

function out=scoreDiv(tcells1,candidates,fr,channel,minDivisionTime)

out=[];


for i=1:length(candidates)
   fluo=[tcells1(candidates(i)).Obj.fluoMean];
   fluo=reshape(fluo,length(tcells1(candidates(i)).Obj(1).fluoMean),[]);
   fluo=fluo(channel,:);
   
   fluo2=[tcells1(candidates(i)).Obj.area];
   
  % size(fluo)size(fluo2)
   fluo=fluo.*fluo2;
   
   
   fluo=smooth(fluo,6);
  % figure, plot(fluo)
   
   fluo=diff(fluo);

   warning off all
   di=min(minDivisionTime,length(fluo)-1);
   [pksmax locmax]=findpeaks(-fluo,'minpeakdistance',di,'minpeakheight',5000); % to be set properly after normalization
   warning on all;
   
  % figure,plot(fluo); hold on; plot(locmax,-pksmax,'Color','r');
  % title(num2str(candidates(i)));
  
  % fr
  % locmax+tcells1(candidates(i)).detectionFrame-1
   
   d=min(abs(locmax+tcells1(candidates(i)).detectionFrame-1-fr));
   
   if d<=2
      out=[out candidates(i)]; 
   end
end


function out=scoreSize2(cells1,candidates)
% identify biggest cells among candidates
out=[];
maxs=-1;
maxind=0;
for i=1:length(cells1)
    % candidates,cells1(i).n
    pix=find(candidates==cells1(i).n);
    if numel(pix)~=0
        % pix
        if cells1(i).area>maxs
            maxs=cells1(i).area;
            maxind=cells1(i).n;
        end
    end
end

out=maxind;

function [out val]=scoreBudNeck(cells1,candidates,targetCell,budneck)

%masks=(zeros(size(displayImage(:,:,1))));

% build mask with budneck label

indbud=[];

for i=1:length(budneck)
    if budneck(i).n~=0
        %        bw_cell = poly2mask(budneck(j).x,budneck(j).y,size(displayImage,1),size(displayImage,2));
        %        masks(bw_cell)=budneck(j).n;
        
        dist=sqrt((budneck(i).ox-targetCell.ox).^2+(budneck(i).oy-targetCell.oy).^2);
        
        if dist<100
           indbud=[indbud i]; 
        end
    end
end

out=[];
val=[];

% identify budneck at the interface between target and candidates

nc=length(cells1);

bud=[];
budval=[];

for i=1:length(candidates)
    
    ind=candidates(i);
    a=[cells1.n];
    pix=find(a==ind);
    
    theta=atan2(cells1(pix).oy-targetCell.oy,cells1(pix).ox-targetCell.ox);
    
    if abs(theta)<pi/4 || abs(theta)>3*pi/4
        xc=[targetCell.ox-3 targetCell.ox-3 cells1(pix).ox+3  cells1(pix).ox+3 targetCell.ox-3];
        yc=[targetCell.oy-3 targetCell.oy+3 cells1(pix).oy+3  cells1(pix).oy-3 targetCell.oy-3];
    else
        xc=[targetCell.ox-3 targetCell.ox+3 cells1(pix).ox+3  cells1(pix).ox-3 targetCell.ox-3];
        yc=[targetCell.oy-3 targetCell.oy-3 cells1(pix).oy+3  cells1(pix).oy+3 targetCell.oy-3];
    end
    
    %ar=polyarea(xc,yc);
    
    inside=[];
    
    for j=indbud
        xb=budneck(j).x;
        yb=budneck(j).y;
       
        frac=double(length(find(inpolygon(xb,yb,xc,yc))))/double(length(xb));
        inside=[inside frac];
    end
    
    [inside ix]=max(inside);
    
    
    if inside>0.1
        bud=[bud ind];
        budval=[budval 30*inside];
    end
    
    %bw_target = poly2mask(xc,yc,size(displayImage,1),size(displayImage,2));
    
    % pix=masks(bw_target);
    
    % if mean(pix)~=0
    %     out=[out ind];
    %     val=[val numel(pix)];
    % end
    
    %inte=find(bw_target & masks);
    
    %if targetCell.n==35
    %ar,theta
    %   figure, imshow(masks,[]); line(xc,yc);
    % end
    
end

out=bud;
val=budval;

%if targetCell.n==35
%out, val
%end


function [problems energy]=computeEnergy(mothers,minDivisionTime,maxDivisionTime,firstframe)

problems=[];

energy=0;

indmothers=[mothers.n];

for i=1:numel(mothers)
    
    
    if numel(mothers(i).budTimes)==0
        continue
    end
    
    % first detect timings issues with daughter cells
    timings=[mothers(i).detectionFrame mothers(i).budTimes];
    delta=timings(2:end)-timings(1:end-1);
    
    pixmin=find(delta<=minDivisionTime);
    
    energy=energy+length(pixmin);
    
    pixmin=[pixmin-1 pixmin];
    if mothers(i).detectionFrame==firstframe
    pixmin=setdiff(pixmin,[0 1]);
    else
    pixmin=setdiff(pixmin,0);    
    end
    
    problems=[problems mothers(i).daughterList(pixmin)];
    
    pixmax=find(delta>=maxDivisionTime);
    
    energy=energy+length(pixmax);
    
    pixmax=[pixmax-1 pixmax];
    pixmax=setdiff(pixmax,0);
    %pixmax=setdiff(pixmax,length(mothers(i).daughterList+1));

    %i,pixmax,mothers(i).daughterList
    %i,pixmin,pixmax,mothers(i).daughterList(pixmin),mothers(i).daughterList(pixmax)
    %pause
    problems=[problems mothers(i).daughterList(pixmax)];
    
    
end

    
    
    
    








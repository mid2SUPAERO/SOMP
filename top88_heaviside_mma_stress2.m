%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
%top88_heaviside(150,50,0.3,1.5,3)
function top88_heaviside_mma_stress2(nelx,nely,rmin,ft,P,Sl,BC,at)
%% MATERIAL PROPERTIES
close all
E0 = 1;
Emin = 1e-9;
nu = 0.3;
penal=3;
beta = 1;
aggregation_settings.ka=P;
aggregation_settings.zp=1;
aggregation_settings.aggregation=at;
% Sl=1;
folder_name=['MMA_Optimization_history_',BC,'nelx_',num2str(nelx),'nely_',num2str(nely),'_R_',num2str(rmin),'_P_',num2str(P),'_ft_',num2str(ft),'_Sl_',num2str(Sl),'_',aggregation_settings.aggregation];
image_prefix=[BC,'nelx_',num2str(nelx),'nely_',num2str(nely),'_R_',num2str(rmin),'_P_',num2str(P),'_ft_',num2str(ft),'_Sl_',num2str(Sl),'_',aggregation_settings.aggregation];
mkdir(folder_name)
Path=[folder_name,'/'];
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
D=E0/(1-nu^2)*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
B=2*nelx/100/2*[-1 0 1 0 1 0 -1 0;0 -1 0 -1 0 1 0 1;-1 -1 -1 1 1 1 1 -1];
DB=D*B;
Cvm=[1 -0.5 0;-0.5 1 0;0 0 3];
Sel=DB'*Cvm*DB;
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
is= reshape(repmat(1:(nelx*nely),8,1),[],1);
js= reshape(edofMat',[],1);
%define the nodal coordinates
[Yy,Xx]=find(nodenrs);
Yy=nely+1-Yy;
Xx=Xx-1;
% Element connectivity
enodeMat=edofMat(:,[2,4,6,8])/2;
% compute the centroid coordinates
xc=mean(Xx(enodeMat'));
yc=mean(Yy(enodeMat'));
centroid_coordinate=[xc(:),yc(:)];
%% DEFINE LOADS AND SUPPORTS
switch BC
    case 'MBB'
        excitation_node=1;excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=[find(Xx==min(Xx));(nelx+1)*(nely+1)];fixed_dir=[ones(nely+1,1);2];
        fixeddofs=2*(fixednodes-1)+fixed_dir;
        emptyelts=[]; fullelts = [];
    case 'Short_Cantilever'
        excitation_node=find((Xx==max(Xx))&(Yy==fix(0.5*min(Yy)+0.5*max(Yy))));excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=repmat(find(Xx==min(Xx)),2,1);fixed_dir=[ones(nely+1,1);2*ones(nely+1,1)];
        fixeddofs=2*(fixednodes-1)+fixed_dir(:);
        emptyelts=[]; fullelts = [];
    case 'L-shape'
        excitation_node=find((Xx==max(Xx))&(Yy==fix(0.5*min(Yy)+0.5*max(Yy))));excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=find(Yy==max(Yy));
        fixednodes=repmat(intersect(fixednodes,find(Xx<=(((max(Xx)+min(Xx))/2)))),2,1);
        fixed_dir=[ones(length(fixednodes)/2,1),2*ones(length(fixednodes)/2,1)];
        fixeddofs=2*(fixednodes-1)+fixed_dir(:);
        nds=0.1;
        fullelts=find(xc>=((((1-nds)*max(Xx)+nds*min(Xx))))&(yc<=((max(Yy)+min(Yy))/2))&(yc>=(((0.5-nds)*max(Yy)+(0.5+nds)*min(Yy)))));
        emptyelts =find(xc>=(((max(Xx)+min(Xx))/2))&(yc>=((max(Yy)+min(Yy))/2)));
    otherwise
        error('BC string should be a valid entry: ''MBB'',''L-Shape'',''Short_Cantilever''')
end
activeelts=setdiff(1:nelx*nely,emptyelts);
U=F;
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
figure(2)
%         map=colormap(gray);
%         map=map(end:-1:1,:);
%         caxis([0 1])
%         patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],'FaceColor','flat','EdgeColor','none');
axis equal; axis off; ;hold on
hold on
if strcmp(BC,'L-shape')
    fill([min(Xx),max(Xx),max(Xx),(min(Xx)+max(Xx))/2,(min(Xx)+max(Xx))/2,min(Xx)],[min(Yy),min(Yy),(min(Yy)+max(Yy))/2,(min(Yy)+max(Yy))/2,max(Yy),max(Yy)],'w')
    fill([(1-nds)*max(Xx)+nds*min(Xx),max(Xx),max(Xx),(1-nds)*max(Xx)+nds*min(Xx)],[(0.5-nds)*max(Yy)+(0.5+nds)*min(Yy),(0.5-nds)*max(Yy)+(0.5+nds)*min(Yy),(min(Yy)+max(Yy))/2,(min(Yy)+max(Yy))/2],'k')
    text(mean(Xx)/2,mean(Yy)/2,'\Omega','FontSize',24,'FontWeight','bold')
else
    fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w')
    text(mean(Xx),mean(Yy),'\Omega','FontSize',24,'FontWeight','bold')
end
scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
scal=10;
quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,-(excitation_direction==2),scal,'r','Linewidth',2)
%     title('density plot')
%         colormap(map)
%         colorbar
hold off
%         axis([min(Xx),max(Xx),min(Yy),max(Yy)])
FEM_enhanceDisplay
drawnow
% print([Path,'design_problem'],'-dpng')
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = ones(nely,nelx);
xPhys = x;
outit = 0;
change = 1;
m = 1;
n = length(xPhys(:));
epsimin = 0.0000001;
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = xPhys(:);
xold1   = xval;
xold2   = xval;
xmin    = zeron;
xmax    = eeen;
low     = xmin;
upp     = xmax;
C       = 1000*eeem;
d       = 0*eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 1000;
kkttol  =0.001;
%
%%%% The iterations start:
kktnorm = kkttol+10;

if ft == 1 || ft == 2
    xPhys = x;
elseif ft == 3
    xTilde = x;
    xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
end
outit = 0;
outitbeta=0;
outitp=0;
outitr=0;
change = 1;
plot_rate=100;
transition=1000;
cvec=zeros(maxoutit,1);
vvec=cvec;
%% START ITERATION
while outit < maxoutit && change>0.001
    outit = outit + 1;
    outitbeta=outitbeta+1;
    outitp=outitp+1;
    outitr=outitr+1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    dK=decomposition(K(freedofs,freedofs));
%     U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    U(freedofs) =dK\F(freedofs);

    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    o=mean(xPhys(:));
    do=ones(size(xPhys))/length(xPhys(:));
    
    mse = reshape(sqrt(abs(sum((U(edofMat)*Sel).*U(edofMat),2))),nely,nelx); %microscopic Von Mises Stress
    rcv=xPhys(:).*(mse(:)/Sl-1);
    rcv(fullelts)=0;
%     rcv(xPhys<=0.015)=0;
    sS=repmat(xPhys(:)'./mse(:)'/Sl,8,1).*(U(edofMat)*Sel)';
    drcv_dU=sparse(is,js,sS,nelx*nely,length(U));
    [c,dGks_drcv]=Aggregation_Pi(rcv,aggregation_settings);
%     dGks_drcv=exp(P*(rcv(:)-max(rcv(:))))/sum(exp(P*(rcv(:)-max(rcv(:)))));
     dG_ksdu=(dGks_drcv'*drcv_dU)';
    
%     sS=reshape(Sel(:)*(xPhys(:)'./mse(:)'/Sl.*exp(P*(rcv(:)'-max(rcv(:)))))/sum(exp(P*(rcv(:)-max(rcv(:))))),64*nelx*nely,1);
%     S0=sparse(iK,jK,sS); S0=(S0+S0')/2;
%     dG_ksdu2=S0*U;
%     norm(dG_ksdu2-dG_ksdu,'inf')
    Lambda=U;
%     Lambda(freedofs) = K(freedofs,freedofs)\dG_ksdu(freedofs);
    Lambda(freedofs) = dK\dG_ksdu(freedofs);

    
    dG_ksdx=(mse(:)/Sl-1).*dGks_drcv;
    % se=mse.*(Emin+reshape(xPhys,nely,nelx).^penal*(E0-Emin))/E0;
%     c=max(rcv(:))+1/P*log(sum(exp(P*(rcv(:)-max(rcv(:))))))-log(length(rcv(:)))/P;%-log(length(rcv(:)))/P
    dc_du=reshape((sum((U(edofMat)*KE).*Lambda(edofMat),2)),nely,nelx);
    dc_du =- penal*(E0-Emin)*xPhys(:).^(penal-1).*dc_du(:);
    dc=dc_du(:)+dG_ksdx;
%     dc(xPhys<=0.015)=0;
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        do(:) = H*(do(:)./Hs);
        dc(:) = H*(dc(:)./Hs);
    elseif ft == 3
        dx = beta*exp(-beta*xTilde)+exp(-beta);
        do(:) = H*(do(:).*dx(:)./Hs);
        dc(:) = H*(dc(:).*dx(:)./Hs);
    end
    dc(emptyelts)=0;
    dc(fullelts)=0;
    do(emptyelts)=0;
    do(fullelts)=0;
    f0val=o;
    fval=c*100;%;1*(-0.38-FAN_Ax_disp)/0.38
    df0dx=do(:);
    dfdx=dc(:)'*100;%;-1*dfandisp(:)'/0.38
    innerit=0;
    outvector1 = [outeriter innerit f0val fval'];
    outvector2 = xval;
    %% MMA code optimization
    [x(:),ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub_conservative(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,C,d);
    xold2 = xold1;
    xold1 = xval;
    xval  = x(:);
    change=norm(xval-xold1,'inf');
    %     xval(x1>=2500&x4<=4500&z1==min(z1))=0;
    if ft == 1
        xPhys(:)=xval;
    elseif ft == 2
        xPhys(:) = (H*xval(:))./Hs;
    elseif ft == 3
        xTilde(:) = (H*xval(:))./Hs;
        xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
    end
    xPhys(emptyelts) = 0;
    xPhys(fullelts) = 1;
     if rmin~=2&&(outitr >= transition|| change <= 0.001)       
        rmin=max(rmin-1,2);
        %% PREPARE FILTER
        iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
        jH = ones(size(iH));
        sH = zeros(size(iH));
        k = 0;
        for i1 = 1:nelx
            for j1 = 1:nely
                e1 = (i1-1)*nely+j1;
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                    end
                end
            end
        end
        H = sparse(iH,jH,sH);
        Hs = sum(H,2);
        outitr=0;
        change=1;
        fprintf('Parameter rmin increased to %g.\n',rmin);  
    end
    if penal<6&&(outitp >= transition || change <= 0.001)
        penal=penal+1;
        
        outitp=0;
        change=1;
        fprintf('Parameter penal increased to %g.\n',penal);
    end
    if ft == 3 && beta < 10 && ...
            (outitbeta >= transition|| change <= 0.001 )%|| change <= 0.01
        beta = beta+1;
        outitbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',outit,full(c), ...
        mean(xPhys(:)),full(change));
    cvec(outit)=c;vvec(outit)=mean(xPhys(:));
    %% PLOT DENSITIES
    if rem(outit,plot_rate)==0
        figure(3)
        subplot(2,1,1)
        plot(1:outit,cvec(1:outit),'bo','MarkerFaceColor','b')
        % title(['Convergence  Compliance =',num2str(c,'%4.3e'),', iter = ', num2str(outit)])
        grid on
        hold on
        scatter(outit,c,'k','fill')
        hold off
        text(outit,c,['C =',num2str(c,'%4.2f'),' at iteration ', num2str(outit)],'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
        % legend(['C =',num2str(c,'%4.3e'),', iter = ', num2str(outit)],'Location','Best')
        xlabel('iter')
        ylabel('C')
        FEM_enhanceDisplay
        % axis([0 outit 0 1e3])
        subplot(2,1,2)
        plot(1:outit,vvec(1:outit)*100,'ro','MarkerFaceColor','r')
        grid on
        hold on
        scatter(outit,mean(xPhys(:))*100,'k','fill')
        hold off
        text(outit,mean(xPhys(:))*100,['V = ',num2str(mean(xPhys(:))*100,'%4.2f'),'% at iteration ', num2str(outit)],'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
        % legend(['V = ',num2str(mean(xPhys(:))*100),', iter = ', num2str(outit)],'Location','Best')
        xlabel('iter')
        ylabel('V [%]')
        FEM_enhanceDisplay
        %         print([Path,image_prefix,'convergence_',num2str(outit-1,'%03d')],'-dpng')
        
        figure(1)
        map=colormap(gray);
        map=map(end:-1:1,:);
        caxis([0 1])
        patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],'FaceColor','flat','EdgeColor','none'); axis equal; axis off; ;hold on
        hold on
        if strcmp(BC,'L-shape')
            fill([min(Xx),max(Xx),max(Xx),(min(Xx)+max(Xx))/2,(min(Xx)+max(Xx))/2,min(Xx)],[min(Yy),min(Yy),(min(Yy)+max(Yy))/2,(min(Yy)+max(Yy))/2,max(Yy),max(Yy)],'w','FaceAlpha',0.)
        else
            fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
        end
%         fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
        scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
        scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
        scal=10;
        quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,-(excitation_direction==2),scal,'r','Linewidth',2)
        %     title('density plot')
        colormap(map)
        colorbar
        hold off
        %         axis([min(Xx),max(Xx),min(Yy),max(Yy)])
        FEM_enhanceDisplay
        drawnow
        %         print([Path,'density_',num2str(outit-1,'%03d')],'-dpng')
    end
    %   colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
figure(1)
map=colormap(gray);
map=map(end:-1:1,:);
caxis([0 1])
patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],'FaceColor','flat','EdgeColor','none'); axis equal; axis off; ;hold on
hold on
if strcmp(BC,'L-shape')
            fill([min(Xx),max(Xx),max(Xx),(min(Xx)+max(Xx))/2,(min(Xx)+max(Xx))/2,min(Xx)],[min(Yy),min(Yy),(min(Yy)+max(Yy))/2,(min(Yy)+max(Yy))/2,max(Yy),max(Yy)],'w','FaceAlpha',0.)
        else
            fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
end
scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
scal=10;
quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,-(excitation_direction==2),scal,'r','Linewidth',2)
%     title('density plot')
colormap(map)
colorbar
drawnow
hold off
% axis([min(Xx),max(Xx),min(Yy),max(Yy)])
FEM_enhanceDisplay
print([Path,image_prefix,'density'],'-dpng')
figure(3)
subplot(2,1,1)
plot(1:outit,cvec(1:outit),'bo','MarkerFaceColor','b')
% title(['Convergence  Compliance =',num2str(c,'%4.3e'),', iter = ', num2str(outit)])
grid on
hold on
scatter(outit,c,'k','fill')
hold off
text(outit,c,['C =',num2str(c,'%4.2f'),' at iteration ', num2str(outit)],'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
% legend(['C =',num2str(c,'%4.3e'),', iter = ', num2str(outit)],'Location','Best')
xlabel('iter')
ylabel('C')
FEM_enhanceDisplay
% axis([0 outit 0 1e3])
subplot(2,1,2)
plot(1:outit,vvec(1:outit)*100,'ro','MarkerFaceColor','r')
grid on
hold on
scatter(outit,mean(xPhys(:))*100,'k','fill')
hold off
text(outit,mean(xPhys(:))*100,['V = ',num2str(mean(xPhys(:))*100,'%4.2f'),'% at iteration ', num2str(outit)],'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
% legend(['V = ',num2str(mean(xPhys(:))*100),', iter = ', num2str(outit)],'Location','Best')
xlabel('iter')
ylabel('V [%]')
FEM_enhanceDisplay
% axis([0 outit 0 Inf])
% title(['Convergence volfrac = ',num2str(mean(xPhys(:))*100),', iter = ', num2str(outit)])
print([Path,image_prefix,'convergence'],'-dpng')
mse = reshape(sqrt(sum((U(edofMat)*Sel).*U(edofMat),2)),nely,nelx); %microscopic Von Mises Stress
mse(fullelts)=0;
rcv=xPhys.*(mse/Sl-1);%relaxed constraint violation
se=mse.*(Emin+xPhys.^penal*(E0-Emin))/E0; %macroscopic stress
cv=se/Sl-1;% constraint violation
figure(4)
map=colormap(jet);
% map=map(end:-1:1,:);
%caxis([0 1])
patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(activeelts(:),[2,4,6,8])/2,'FaceVertexCData',full(se(activeelts(:))),'FaceColor','flat','EdgeColor','none'); axis equal; axis off; ;hold on
hold on
if strcmp(BC,'L-shape')
            fill([min(Xx),max(Xx),max(Xx),(min(Xx)+max(Xx))/2,(min(Xx)+max(Xx))/2,min(Xx)],[min(Yy),min(Yy),(min(Yy)+max(Yy))/2,(min(Yy)+max(Yy))/2,max(Yy),max(Yy)],'w','FaceAlpha',0.)
        else
            fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
end
scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
scal=10;
quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,-(excitation_direction==2),scal,'r','Linewidth',2)
%     title('density plot')
colormap(map)
colorbar
drawnow
hold off
% axis([min(Xx),max(Xx),min(Yy),max(Yy)])
FEM_enhanceDisplay
print([Path,image_prefix,'VM_stress'],'-dpng')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


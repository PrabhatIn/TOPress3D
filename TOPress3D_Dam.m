function TOPress3D_Dam(nelx,nely,nelz,volfrac,penal,rmin,eta,beta,lst,maxit)
%% ___PART 1.____________________________MATERIAL AND FLOW PARAMETERS
 E1 = 1;                                                                    % Youngs' Modulus of solid
 Emin = E1*1e-5;                                                            % Youngs's Modulus of void
 nu = 0.30;                                                                 % Poisson's ratio
[Kv,epsf,r,Dels] = deal(1,1e-7,0.1,2);                                     % Flow parameters
[Ds, Kvs]= deal((log(r)/Dels)^2*epsf,Kv*(1 - epsf));                       % Flow parameters 
%% ____PART 2._______________FINITE ELEMENT ANALYSIS PREPARATION and NON-DESIGN REGIONS
[ndx,ndy,ndz] = deal(nelx+1,nely+1,nelz+1);                                % Node number in x, y and z
[nel,nno] = deal(nelx*nely*nelz, ndx*ndy*ndz);                             % Number of elements and nodes
nodenrs = int32(reshape( 1 : nno,ndy, ndz, ndx ) );   
edofVec = reshape( 3 * nodenrs( 1 : nely, 1 : nelz, 1 : nelx ) + 1, nel, 1);
Udofs = edofVec + int32( [0,1,2,3*ndy*ndz+[0,1,2,-3,-2,-1],-3,-2,-1,3*ndy+...
    [0,1,2],3*ndy*(ndz+1)+[0,1,2,-3,-2,-1],3*ndy+[-3,-2,-1]]);             % Displacement DOFs matrix
[Pdofs,allPdofs, allUdofs] = deal(Udofs(:,3:3:end)/3,1:nno,1:3*nno);       % Pressure DOFs matrix, all DOFs and displacement DOFs
Tface = (ndy*nelz+1:ndz*ndy:nno-nely)+ (0:nely)'; 
BTface=(1:ndz*ndy:nno-nely)+ (0:nely)';
[Tface, BTface] = deal(Tface(:)', BTface(:)');
[Lface,Rface] = deal((1:ndy*ndz), (ndy*ndz*nelx+1:nno));                   % Determinig left and right faces
[Bface,Fface] = deal((ndy:ndy:nno),(1:ndy:nno-nely));                      % Determining front and 
[skI, skII, spI, spII] = deal( [ ] );                                       
for j = 1 : 24
    skI = cat( 2, skI, j : 24 );
    skII = cat( 2, skII, repmat( j, 1, 24 - j + 1 ) );
end
[iK, jK ] = deal( Udofs( :,  skI )', Udofs( :, skII )' );
Iar = sort( [ iK( : ), jK( : ) ], 2, 'descend' ); clear iK jK              % Reduced assembly indexing
Ke = 1/(1+nu)/(2*nu-1)/144 *( [ -32;-6;-6;8;6;6;10;6;3;-4;-6;-3;-4;-3;-6;10;...
    3;6;8;3;3;4;-3;-3; -32;-6;-6;-4;-3;6;10;3;6;8;6;-3;-4;-6;-3;4;-3;3;8;3;...
    3;10;6;-32;-6;-3;-4;-3;-3;4;-3;-6;-4;6;6;8;6;3;10;3;3;8;3;6;10;-32;6;6;...
    -4;6;3;10;-6;-3;10;-3;-6;-4;3;6;4;3;3;8;-3;-3;-32;-6;-6;8;6;-6;10;3;3;4;...
    -3;3;-4;-6;-3;10;6;-3;8;3;-32;3;-6;-4;3;-3;4;-6;3;10;-6;6;8;-3;6;10;-3;...
    3;8;-32;-6;6;8;6;-6;8;3;-3;4;-3;3;-4;-3;6;10;3;-6;-32;6;-6;-4;3;3;8;-3;...
    3;10;-6;-3;-4;6;-3;4;3;-32;6;3;-4;-3;-3;8;-3;-6;10;-6;-6;8;-6;-3;10;-32;...
    6;-6;4;3;-3;8;-3;3;10;-3;6;-4;3;-6;-32;6;-3;10;-6;-3;8;-3;3;4;3;3;-4;6;...
    -32;3;-6;10;3;-3;8;6;-3;10;6;-6;8;-32;-6;6;8;6;-6;10;6;-3;-4;-6;3;-32;6;...
    -6;-4;3;6;10;-3;6;8;-6;-32;6;3;-4;3;3;4;3;6;-4;-32;6;-6;-4;6;-3;10;-6;3;...
    -32;6;-6;8;-6;-6;10;-3;-32;-3;6;-4;-3;3;4;-32;-6;-6;8;6;6;-32;-6;-6;-4;...
    -3;-32;-6;-3;-4;-32;6;6;-32;-6;-32]+nu*[ 48;0;0;0;-24;-24;-12;0;-12;0;...
    24;0;0;0;24;-12;-12;0;-12;0;0;-12;12;12;48;0;24;0;0;0;-12;-12;-24;0;-24;...
    0;0;24;12;-12;12;0;-12;0;-12;-12;0;48;24;0;0;12;12;-12;0;24;0;-24;-24;0;...
    0;-12;-12;0;0;-12;-12;0;-12;48;0;0;0;-24;0;-12;0;12;-12;12;0;0;0;-24;...
    -12;-12;-12;-12;0;0;48;0;24;0;-24;0;-12;-12;-12;-12;12;0;0;24;12;-12;0;...
    0;-12;0;48;0;24;0;-12;12;-12;0;-12;-12;24;-24;0;12;0;-12;0;0;-12;48;0;0;...
    0;-24;24;-12;0;0;-12;12;-12;0;0;-24;-12;-12;0;48;0;24;0;0;0;-12;0;-12;...
    -12;0;0;0;-24;12;-12;-12;48;-24;0;0;0;0;-12;12;0;-12;24;24;0;0;12;-12;...
    48;0;0;-12;-12;12;-12;0;0;-12;12;0;0;0;24;48;0;12;-12;0;0;-12;0;-12;-12;...
    -12;0;0;-24;48;-12;0;-12;0;0;-12;0;12;-12;-24;24;0;48;0;0;0;-24;24;-12;...
    0;12;0;24;0;48;0;24;0;0;0;-12;12;-24;0;24;48;-24;0;0;-12;-12;-12;0;-24;...
    0;48;0;0;0;-24;0;-12;0;-12;48;0;24;0;24;0;-12;12;48;0;-24;0;12;-12;-12;...
    48;0;0;0;-24;-24;48;0;24;0;0;48;24;0;0;48;0;0;48;0;48 ] );             % Elemental stiffness matrix  
Ke0( tril( ones( 24 ) ) == 1 ) = Ke';
Ke0 = reshape( Ke0, 24, 24 ); 
Ke0 = Ke0 + Ke0' - diag( diag( Ke0 ) );                                    % Extracting full form
for j = 1 : 8
    spI = cat( 2, spI, j : 8 );
    spII = cat( 2, spII, repmat( j, 1, 8-j + 1 ) );
end
[iP, jP ] = deal( Pdofs( :,  spI )', Pdofs( :, spII )' );
IarP = sort( [ iP( : ), jP( : ) ], 2, 'descend' ); clear iP jP   
Kpl = [4;0;-1;0;0;-1;-1;-1;4;0;-1;-1;0;-1;-1;4;0;-1;-1;0;-1;4;-1;-1;-1; ...
            0;4;0;-1;0;4;0;-1;4;0;4]/12;                                      % Flow matrix due to Darcy
KDpl = [8;4;2;4;4;2;1;2;8;4;2;2;4;2;1;8;4;1;2;4;2;8;2;1;2;4;8;4;2;4;8;...
             4;2;8;4;8]/216;                                                   % Flow matrix due to drainage term
[Kp(tril( ones(8))== 1), KDp(tril(ones(8))==1)]= deal(Kpl',KDpl);
[Kp,KDp]= deal(reshape( Kp, 8, 8 ),reshape( KDp, 8, 8 ));
Kp = Kp + Kp' - diag( diag( Kp ) ); KDp = KDp + KDp' - diag( diag( KDp ) );% Extracting full form
[lKe, lKpl] = deal(length(Ke), length(Kpl));
Te = [-4;-4;-4;-4;-2;-2;-2;	-2;	-1;	-2;-4;-2;-2;-2;-4;-2;-1;-2;-1;-1;-1;-1;...
   -2;-2;4;-2;-2;4;-4;-4;2;-4;-2;2;-2;-1;2;-1;-2;2;-2;-4;1;-2;-2;1;-1;-1;2;...
   2;-1;2;4;-2;4;4;-4;4;2;-2;1;1;-1;1;2;-2;2;2;-4;2;1;-2;-2;4;-2;-2;2;-1;-4;...
   2;-2;-4;4;-4;-1;2;-2;-1;1;-1;-2;1;-2;-2;2;-4;-2;-2;4;-2;-1;2;-1;-1;1;-1;...
   -2;2;-4;-4;4;-4;-2;2;-2;-2;1;-2;-4;2;2;-1;2;2;-2;4;1;-2;2;1;-1;1;4;-2;2;...
   4;-4;4;2;-4;2;2;-2;1;1;1;1;1;2;2;2;2;4;2;1;2;2;2;1;2;4;2;4;4;4;4;2;2;-1;...
   2;2;-1;1;1;-2;1;2;-2;2;4;-2;4;2;-2;2;1;-4;2;2;-4;4;4]/72;
iT = reshape(kron(Udofs,int32(ones(8,1)))',192*nel,1);
jT = reshape(kron(Pdofs,int32(ones(1,24)))',192*nel,1);
Ts = reshape(Te(:)*ones(1,nel), 192*nel, 1);                               % Elemental transformation matrix
TG = fsparse(iT, jT, Ts); clear Te iT jT Ts                                                % Global transformation matrix
IFprj=@(xv,etaf,betaf)((tanh(betaf*etaf) + tanh(betaf*(xv-etaf)))/...      % Projection function
    (tanh(betaf*etaf) + tanh(betaf*(1 - etaf))));
dIFprj=@(xv,etaf,betaf) betaf*(1-tanh(betaf*(xv-etaf)).^2)...
    /(tanh(betaf*etaf)+tanh(betaf*(1-etaf)));                              % Derivative of the projection function
[NDS, NDV ] = deal( [], [] );
act = setdiff((1 : nel)', union( NDS, NDV ));
%% ____PART 3.______PRESSURE & STRUCTURE B.C's, LOADs
[PF, Pin] =deal(0.00001*ones(nno,1),1);                                    % Pressure-field preparation
PF(Fface) = 0; PF(Bface) = Pin;         % Applying pressure load
fixedPdofs = allPdofs(PF~=0.00001);                                        % Given P-dofs
freePdofs  = setdiff(allPdofs,fixedPdofs);                                 % Free P-dofs
pfixeddofsv = [fixedPdofs' PF(fixedPdofs)];                                % p-fixed and its value
fixnn = unique([BTface,Rface]);
fixedUdofs = [3*fixnn-2  3*fixnn-1  3*fixnn  3*Lface-2];                         % Fixed displ.
freeUdofs = setdiff(allUdofs,fixedUdofs);                                  % Free dofs for displ.
[U, lam1] = deal(zeros(3*nno,1),zeros(nno,1));                             % lam1:Lagrange mult.
%% ___PART 4._________________________________________FILTER PREPARATION
[dy,dz,dx]=meshgrid(-ceil(rmin)+1:ceil(rmin)-1,...
    -ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1 );
h = max( 0, rmin - sqrt( dx.^2 + dy.^2 + dz.^2 ) );                        % Conv. kernel                
Hs = imfilter( ones( nely, nelz, nelx ), h);                               % Matrix of weights (filter)  
%% ___PART 5.__________________________MMA OPTIMIZATION PREPARATION & INITIALIZATION
[x,dVol0] = deal(zeros(nel,1),ones(nel,1)/(nel*volfrac));                  % Design var. and vol. cont. der.
 x(act) = (volfrac*(nel-length(NDV))-length(NDS) )/length(act); x(NDS) = 1;% Updating
[nMMA,mMMA,xphys,xMMA,mvLt] = deal(length(act),1,x,x(act),0.1);            % Different variables
[xminvec,xmaxvec] = deal(zeros(nMMA,1),ones(nMMA,1));                      % Min. & Max vector for MMA
[low, upp, xold1,xold2] = deal(xminvec,xmaxvec,xMMA,xMMA);                                        % Low and Upp limits MMA
[cMMA,dMMA, a0, aMMA] = deal(1000*ones(mMMA,1),zeros(mMMA,1),1,zeros(mMMA,1));                                       
dVol = imfilter(reshape(dVol0, nely, nelz, nelx)./Hs,h);                   % Filtered volume sensitivity
[loop, change] =deal(0,1);   
%% ____PART 6._____________________________________MMA OPTIMIZATION LOOP
while(loop<maxit && change>0.0001)
 loop = loop + 1;                                                          % Updating the opt. iteration
 %___PART 6.1__________SOLVING FLOW BALANCE EQUATION
Kc = Kv*(1-(1-epsf)*IFprj(xphys,eta,beta));                                % Flow coefficient
Dc = Ds*IFprj(xphys,eta,beta);                                             % Drainage coefficient
Ae = reshape(Kpl(:)*Kc' + KDpl(:)*Dc',lKpl*nel,1);                         % Elemental flow matrix
AG = fsparse(IarP(:,1),IarP(:,2),Ae,[nno, nno]) ;                          % Global flow matrix
Aff = AG(freePdofs,freePdofs);                                             % AG for free pressure dofs
AG = AG + AG' - diag( diag( AG ) ) ;                                       % Determining full AG
PF(freePdofs) = decomposition(Aff,'ldl','lower')\(-AG(freePdofs,fixedPdofs)*pfixeddofsv(:,2)); % Solving for pressure
PF(pfixeddofsv(:,1)) = pfixeddofsv(:,2);                                   % Final P-field
%__PART 6.2_DETERMINING CONSISTENT NODAL LOADS and GLOBAL Disp. Vector
F = -TG*PF;                                                                % Dertmining nodal forces
E = Emin + xphys.^penal*(E1 - Emin);                                       % Material interpolation
Ks = reshape(Ke(:)*E',lKe*nel,1);                                          % Elemental stiffness matrix
KG = fsparse( Iar( :, 1 ), Iar( :, 2 ), Ks, [3*nno, 3*nno ] );             % Global stiffnes matrix
L = chol( KG( freeUdofs, freeUdofs ), 'lower' );
U(freeUdofs ) = L' \ ( L \ F( freeUdofs ) );   
%__PART 6.3__OBJECTIVE, CONSTRAINT and THEIR SENSITIVITIES COMPUTATION
obj = U'*F;                                                                % Determining objective
lam1(freePdofs) = (2*U(freeUdofs)'*TG(freeUdofs,freePdofs))/Aff;           % Lagrange mult.
objsT1 = -(E1 - Emin)*penal*xphys.^(penal - 1).*sum(([U(Udofs)]*Ke0).*[U(Udofs)],2);
dC1k = -dIFprj(xphys,eta,beta).* sum((lam1(Pdofs)*(Kvs*Kp)) .* PF(Pdofs),2);
dC1d =  dIFprj(xphys,eta,beta).* sum((lam1(Pdofs)*(Ds*KDp)) .* PF(Pdofs),2);
objsT2 = dC1k + dC1d; 
objsens = (objsT1 + lst*objsT2);                                           % Final sensitivities
Vol = sum(xphys)/(nel*volfrac)-1;                                          % Volume fraction
if(loop ==1), normf = 1000/(obj);end
objsens = imfilter(reshape(objsens*normf, nely, nelz, nelx)./Hs,h);        % Obj. sens.
%___PART 6.4______________________SETTING and CALLING MMA OPTIMIZATION
xval = xMMA;
[xminvec, xmaxvec]= deal(max(0, xval - mvLt),min(1, xval + mvLt));
[xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(mMMA,nMMA,loop,xval,xminvec,xmaxvec,...
    xold1,xold2,obj*normf,objsens(act),Vol,dVol(act)',low,upp,a0,aMMA,cMMA,dMMA);                               % calling MMA
[xold2,xold1,xnew] = deal(xold1, xval, xmma);                              % Updating
change = max(abs(xnew-xMMA)); xMMA = xnew;                                 % Updating solutions                                              % Calculating chan
xphys(act) = xnew; 
xphys = imfilter(reshape(xphys, nely, nelz, nelx),h)./Hs;% Filering the Phy. vector
xphys = xphys(:); xphys(NDS)= 1; xphys(NDV)= 0;                              % Updating xphys for active/passive region
%___PART 6.5_____________________________PRINTING and PLOTTING RESULTS
fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,obj*normf,mean(xphys),change);
cla;isovals = zeros(nelx*2,nely,nelz);
isovals(nelx+1:2*nelx,1:nely,1:nelz) = shiftdim(reshape(xphys,nely,nelz,nelx),2);
isovals(1:nelx,1:nely,1:nelz) = isovals(2*nelx:-1:nelx+1,1:nely,1:nelz);
isovals = smooth3( isovals, 'box', 1 );
patch(isosurface(isovals, 0.30),'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
patch(isocaps(isovals, 0.30),'FaceColor','r','EdgeColor','none');
view( 3 ); axis equal tight off; drawnow, camlight; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%    The code, TOPress3D,  is provided for  pedagogical purposes. A  detailed        %
	%    description is presented in the paper:"TOPress3D: 3D topology optimization      %       
	%    with design-dependent pressure loads in MATLAB"   Optimization and Engineering, % 
	%     2024.                                                                          %
	%    One can download the code and its extensions for the different problems         %
	%    from the online supplementary material and also from:                           %
	%                                      https://github.com/PrabhatIn/TOPress3D        %
	%    Please send your comment to: pkumar@mae.iith.ac.in                              %
	%    One may also refer to the following four papers:                                % 
	%                                                                                    %
	%    1. Kumar P, Frouws JS, Langelaar M (2020) Topology optimization of fluidic      %
	%    pressure-loaded structures and compliant mechanisms using the Darcy method.     %
	%    Structural and Multidisciplinary Optimization 61(4):1637-1655                   %
	%    2. Kumar P, Langelaar M (2021) On topology optimization of design-dependent     % 
	%    pressure-loaded three-dimensional structures and compliant mechanisms.          %
	%    International Journal for Numerical Methods in Engineering 122(9):2205-2220     %
	%    3. Kumar P. (2023) TOPress: a MATLAB implementation for topology optimization   %
	%    of structures subjected to design-dependent pressure loads.                     % 
	%    Structural and Multidisciplinary Optimization 66(4):97                          %
	%    4. Kumar P. (2024) SoRoTop: a hitchhiker's guide to topology optimization       %
	%    MATLAB code for design-dependent pneumatic-driven soft robots                   % 
	%    Optimization and Engineering, 25:2473-2507(2024)                                %
	%                                                                                    %
	%    Disclaimer:                                                                     %
	%    The author does not guarantee that the code is free from erros but reserves     %
	%    all rights. Further, the author shall not be liable in any event caused by      % 
	%    use of the above mentioned code and its extensions                              %
	%                                                                                    %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

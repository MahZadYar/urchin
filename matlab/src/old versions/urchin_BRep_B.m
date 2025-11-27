function mesh = urchin_BRep(varargin)
% URCHIN_BREP  Structured nano-urchin (Option B explicit spikes)
% Eliminates alphaShape bridging by using ring-based spike construction.
% See header in earlier discussion for parameter details.

p=inputParser; p.KeepUnmatched=false;
addParameter(p,'cr',30,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'sl',15,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'ns',50,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'st',5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'sc',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'filletRatio',0.3,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
addParameter(p,'density',8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'distMethod','uniform');
addParameter(p,'sf',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'flucMethod','uniform');
addParameter(p,'nConeRings',[],@(x)isempty(x)||(isscalar(x)&&x>=4));
addParameter(p,'nTipRings',[],@(x)isempty(x)||(isscalar(x)&&x>=2));
addParameter(p,'nFilletRings',[],@(x)isempty(x)||(isscalar(x)&&x>=1));
parse(p,varargin{:}); S=p.Results;

cr=S.cr; sl=S.sl; ns=S.ns; st=S.st; sc=S.sc; fr=S.filletRatio; dens=S.density;
r_tip=st/2; if ns>1, hcap=2*cr/ns; else, hcap=cr; end
r_base_max=sqrt(max(0,2*cr*hcap-hcap^2)); r_base=max(r_tip*1.05, sc*r_base_max);
fillet_r=fr*r_base; if fillet_r<1e-9||fr<=0, fillet_r=0; end; body_len=sl;

switch lower(S.distMethod)
    case 'uniform'; dirs=fibonacci_sphere(ns,1);
    otherwise; phi=2*pi*rand(ns,1); cz=2*rand(ns,1)-1; th=acos(cz); dirs=[sin(th).*cos(phi),sin(th).*sin(phi),cz];
end
if S.sf>0
    switch lower(S.flucMethod)
        case 'uniform'
            if exist('sobolset','file'), v=net(sobolset(1),ns)-0.5; else, v=rand(ns,1)-0.5; end
        otherwise
            rng(ns,'twister'); v=rand(ns,1)-0.5;
    end
    dL=2*v*S.sf*sl;
else
    dL=zeros(ns,1);
end

V=zeros(0,3); F=zeros(0,3); fprintf('Generating %d spikes...\n',ns);
for i=1:ns
 o=dirs(i,:); [u,v]=local_frame(o); Ls=max(0.2*sl, body_len+dL(i));
 nSeg=max(16,ceil(2*pi*r_base*sqrt(dens/10))); th=linspace(0,2*pi,nSeg+1); th(end)=[];
 nCone=max(6,ceil(10*Ls/sl)); if ~isempty(S.nConeRings), nCone=S.nConeRings; end
 nTip=4; if ~isempty(S.nTipRings), nTip=S.nTipRings; end
 nFil= fillet_r>0 * max(2,ceil(4*fr)); if ~isempty(S.nFilletRings)&&fillet_r>0, nFil=S.nFilletRings; end
 rings=[];
 if fillet_r>0 && nFil>0
  for k=0:nFil
   t=k/nFil; r_here=r_base - fillet_r*(1-cos(t*pi/2)); z_pen=fillet_r*sin(t*pi/2); c=(cr - z_pen)*o;
   pts=circle_points(c,u,v,r_here,th); [V,idx]=append_vertices(V,pts); rings=[rings; idx']; %#ok<AGROW>
  end
 else
  c=cr*o; pts=circle_points(c,u,v,r_base,th); [V,idx]=append_vertices(V,pts); rings=[rings; idx'];
 end
 for k=2:nCone
  t=(k-1)/(nCone-1); r_here=r_base*(1-t)+r_tip*t; c=(cr+t*Ls)*o; pts=circle_points(c,u,v,r_here,th);
  [V,idx]=append_vertices(V,pts); rings=[rings; idx']; %#ok<AGROW>
 end
 for k=1:nTip
  t=k/(nTip+1); ang=(pi/2)*(1-t); rR=r_tip*sin(ang); z_off=r_tip*(1-cos(ang)); c=(cr+Ls+z_off)*o; pts=circle_points(c,u,v,rR,th);
  [V,idx]=append_vertices(V,pts); rings=[rings; idx']; %#ok<AGROW>
 end
 apex=(cr+Ls+r_tip)*o; [V,ap]=append_vertices(V,apex);
 for r=1:size(rings,1)-1, F=[F; connect_rings(rings(r,:),rings(r+1,:))]; end %#ok<AGROW>
 F=[F; connect_ring_to_apex(rings(end,:),ap)]; %#ok<AGROW>
end

fprintf('Meshing core...\n'); nCore=max(400,ceil(dens*4*pi*cr^2/80)); core=fibonacci_sphere(nCore,cr); keep=true(size(core,1),1);
excl=asin(min(0.99,(r_base+fillet_r)/cr)); for i=1:ns, o=dirs(i,:); ang=acos(max(-1,min(1,(core*o')/cr))); keep=keep & (ang>excl); end
core=core(keep,:); try, shp=alphaShape(core(:,1),core(:,2),core(:,3)); a=criticalAlpha(shp,'one-region'); if ~isempty(a), shp.Alpha=a(1)*1.05; end
 [Fc,Vc]=boundaryFacets(shp); b=size(V,1)+1; V=[V;Vc]; F=[F; Fc+b-1]; catch, warning('core alphaShape failed'); end

[Vw,~,ic]=uniquetol(V,1e-7,'ByRows',true); Fw=ic(F); mesh.Vertices=Vw; mesh.Faces=Fw; fprintf('Done V=%d F=%d\n',size(Vw,1),size(Fw,1));
if nargout==0, figure('Color','w'); patch('Faces',Fw,'Vertices',Vw,'FaceColor',[0.9 0.7 0.2],'EdgeColor','none'); axis equal off; view(3); camlight headlight; lighting gouraud; end
end

%% Helpers
function pts=circle_points(c,u,v,r,th); ca=cos(th); sa=sin(th); pts=c + r*(ca'.*u + sa'.*v); end
function [V,idx]=append_vertices(V,P); s=size(V,1)+1; V=[V;P]; idx=(s:size(V,1)).'; end
function F=connect_rings(a,b); n=numel(a); F=zeros(2*n,3); for k=1:n, k2=mod(k,n)+1; F(2*k-1,:)=[a(k),b(k),b(k2)]; F(2*k,:)=[a(k),b(k2),a(k2)]; end; end
function F=connect_ring_to_apex(ring,ap); n=numel(ring); F=zeros(n,3); for k=1:n, k2=mod(k,n)+1; F(k,:)=[ring(k), ring(k2), ap]; end; end
function [u,v]=local_frame(n); n=n(:)/norm(n); if abs(n(1))>0.9, a=[0;1;0]; else, a=[1;0;0]; end; u=cross(n,a); u=u/norm(u); v=cross(n,u); v=v/norm(v); u=u(:)'; v=v(:)'; end
function pts=fibonacci_sphere(N,R); if nargin<2, R=1; end; pts=zeros(N,3); ga=pi*(3-sqrt(5)); for i=0:N-1, y=1-2*(i+0.5)/N; rxy=sqrt(max(0,1-y*y)); th=mod(i*ga,2*pi); pts(i+1,:)=R*[cos(th)*rxy,y,sin(th)*rxy]; end; end

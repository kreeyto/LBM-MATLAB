%% lattice boltzmann multicompoent single distribution
clear all
%% thread-safe LB

nlinks=9;
tau=0.505;
cssq=1.0/3.0;
visc_LB=cssq*(tau-0.5);
visc1=visc_LB;
tau=0.505;
visc_LB=cssq*(tau-0.5);
visc2=visc_LB;
taud=1;
omega=1/tau;

one_ov_nu=1.0/visc_LB;
sharp_c=0.15*3; %*3;%0.12; %D.*1.6; %D./0.6; $sharepning of the interface, antidiffusion term
sigma=0.1;

nx=150;
ny=150;
nsteps=10000;
stamp=10;

ggx=0*10^-6;
ggy=0;

f=zeros(nx,ny,9);
g=zeros(nx,ny,5);
rho=0.*ones(nx,ny);
u=zeros(nx,ny);
v=zeros(nx,ny);
pxx=ones(nx,ny);
pyy=ones(nx,ny);
pxy=ones(nx,ny);

fi=ones(nx,ny).*0;     % current phase field sol
normx=zeros(nx,ny);
normy=zeros(nx,ny);
curvature=zeros(nx,ny);
indicator=zeros(nx,ny);
bool_ind=zeros(nx,ny);
ffx=zeros(nx,ny);
ffy=zeros(nx,ny);
grad_rho_x=zeros(nx,ny);
grad_rho_y=zeros(nx,ny);
mod_grad=zeros(nx,ny);

p=[4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
p_g=[2/6,1/6,1/6,1/6,1/6];
opp=[1, 4, 5, 2, 3, 8, 9, 6, 7];
ex=[0,1,0,-1,0,1,-1,-1,1];
ey=[0,0,1,0,-1,1,1,-1,-1];
fneq=zeros(9,1);
isfluid=zeros(nx,ny);
isfluid(2:nx-1,2:ny-1)=1;

%rhocoarse(nx/2-2:nx/2+2,ny/2-2:ny/2+2)=0;
rho(:,:)=1;
% for ii=1:nx
%     for jj=1:ny
%         if(ii<(nx/2))
%             fi(ii,jj)=1;
%         else
%             fi(ii,jj)=0;
%         end
%     end
% end
% for ii=2:(nx-1)
%     for jj=2:ny-1
%         if((ii-(nx/2))^2+(jj-(ny/2))^2<=30^2)
%             fi(ii,jj)=1;
%         end
%     end
% end
for ii=2:(nx-1)
    for jj=2:ny-1
        
            Ri=sqrt((ii-(nx/2))^2/2.^2+(jj-(ny/2))^2);
            fi(ii,jj) = 0.5 + 0.5 * tanh(2*(20-Ri)/3);
    end
end

f(:,:,1)=p(1)*rho(:,:);
f(:,:,2)=p(2)*rho(:,:);
f(:,:,3)=p(2)*rho(:,:);
f(:,:,4)=p(2)*rho(:,:);
f(:,:,5)=p(2)*rho(:,:);
f(:,:,6)=p(6)*rho(:,:);
f(:,:,7)=p(6)*rho(:,:);
f(:,:,8)=p(6)*rho(:,:);
f(:,:,9)=p(6)*rho(:,:);
g(:,:,1)=p_g(1)*fi(:,:);
g(:,:,2)=p_g(2)*fi(:,:);
g(:,:,3)=p_g(2)*fi(:,:);
g(:,:,4)=p_g(2)*fi(:,:);
g(:,:,5)=p_g(2)*fi(:,:);

for tt=1:nsteps
    %% phase field calc
    for ii=1:nx
        for jj=1:ny
            if(isfluid(ii,jj)==1)
                fi(ii,jj)=sum(g(ii,jj,:),3);
            end
        end 
    end
    for ii=1:nx
        for jj=1:ny
            if(isfluid(ii,jj)==1)
    %% normal calculation and arrays
                grad_fix=0;
                grad_fiy=0;
                for kk=1:9  
                      grad_fix=grad_fix+ 3*p(kk).*ex(kk).*((fi(ii+ex(kk),jj+ey(kk))));
                      grad_fiy=grad_fiy+ 3*p(kk).*ey(kk).*((fi(ii+ex(kk),jj+ey(kk))));
                end
                mod_grad(ii,jj)=sqrt(grad_fix.^2+grad_fiy.^2);
                normx(ii,jj)=grad_fix./(mod_grad(ii,jj)+10^-9);
                normy(ii,jj)=grad_fiy./(mod_grad(ii,jj)+10^-9);
                indicator(ii,jj)=sqrt(grad_fix.^2+grad_fiy.^2);
                            
            end
        end
    end
    %% curvature
    for ii=1:nx
        for jj=1:ny
            if(isfluid(ii,jj)==1)
                 %% surface tension force
                 curvature(ii,jj)=0;
                     for kk=1:9
                        curvature(ii,jj)=curvature(ii,jj) - 3.*p(kk).*(ex(kk).*(normx(ii+ex(kk),jj+ey(kk))) + ...
                                                                       ey(kk).*(normy(ii+ex(kk),jj+ey(kk))));
                     end 
                 ffx(ii,jj)=sigma.*curvature(ii,jj).*normx(ii,jj).*indicator(ii,jj) ;%*bool_ind(ii,jj);
                 ffy(ii,jj)=sigma.*curvature(ii,jj).*normy(ii,jj).*indicator(ii,jj) ;%*bool_ind(ii,jj);
            end
        end
    end
    %% momenti
    for ii=1:nx
        for jj=1:ny
            if(isfluid(ii,jj)==1)
                
                u(ii,jj)=(f(ii,jj,2)-f(ii,jj,4)+f(ii,jj,6)-f(ii,jj,7)-f(ii,jj,8)+f(ii,jj,9))./rho(ii,jj) + ffx(ii,jj)*0.5./rho(ii,jj) ;
                v(ii,jj)=(f(ii,jj,3)-f(ii,jj,5)+f(ii,jj,6)+f(ii,jj,7)-f(ii,jj,8)-f(ii,jj,9))./rho(ii,jj) + ffy(ii,jj)*0.5./rho(ii,jj) ;   
                uu=0.5*(u(ii,jj).^2+v(ii,jj).^2)/cssq;
                rho(ii,jj)=sum(f(ii,jj,:),3) ; 
                for kk=1:9
                    udotc=(u(ii,jj)*ex(kk) + v(ii,jj)*ey(kk))/cssq;
                    
                    HeF=(p(kk)*(rho(ii,jj) + rho(ii,jj).*(udotc + 0.5.*udotc.^2 - uu))).*((ex(kk)-u(ii,jj)).*ffx(ii,jj) + (ey(kk)-v(ii,jj)).*ffy(ii,jj))./(rho(ii,jj).*cssq);
                    
                    feq=p(kk)*(rho(ii,jj) + rho(ii,jj).*(udotc + 0.5.*udotc.^2 - uu)) - 0.5.*(HeF);
                    fneq(kk)=f(ii,jj,kk)- feq;
                end
                pxx(ii,jj)= fneq(2) + fneq(4) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
                pyy(ii,jj)= fneq(3) + fneq(5) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
                pxy(ii,jj)= fneq(6) - fneq(7) + fneq(8) - fneq(9); 
            end
        end 
    end
    %% collision

    for ii=1:nx
        for jj=1:ny
            if(isfluid(ii,jj)==1)
                 uu=0.5*(u(ii,jj).^2+v(ii,jj).^2)/cssq;
                 %visc_loc=fi(ii,jj).*visc1 + (1-fi(ii,jj)).*visc2;
                 %omega=1./(visc_loc*3+0.5);
                 for kk=1:9                                          
                     udotc=(u(ii,jj)*ex(kk) + v(ii,jj)*ey(kk))/cssq;
                     feq=p(kk)*(rho(ii,jj) + rho(ii,jj).*(udotc + 0.5.*udotc.^2 - uu));
                     

                     HeF=0.5.*(p(kk)*(rho(ii,jj) + rho(ii,jj).*(udotc + 0.5.*udotc.^2 - uu))).*((ex(kk)-u(ii,jj)).*ffx(ii,jj) + (ey(kk)-v(ii,jj)).*ffy(ii,jj))./(rho(ii,jj).*cssq);
                    
                     fneq=(ex(kk).*ex(kk)-cssq)*pxx(ii,jj)+(ey(kk).*ey(kk)-cssq)*pyy(ii,jj) ...
                         + 2*ex(kk).*ey(kk).*pxy(ii,jj);

                     f(ii+ex(kk),jj+ey(kk),kk)=feq + (1.0-omega)*(p(kk)/(2*cssq^2))*fneq  + HeF  ;%+guoGrad
                     
                 end

                 for kk=1:5
                     udotc=(u(ii,jj)*ex(kk) + v(ii,jj)*ey(kk))/cssq;
                     feq=p_g(kk).*fi(ii,jj).*(1 + udotc);
                     
                     Hi=sharp_c.*fi(ii,jj).*(1-fi(ii,jj)).*(ex(kk).*normx(ii,jj)+ ey(kk).*normy(ii,jj));%.*bool_ind(ii,jj);
                     %Nci=NCI
                     
                     g(ii,jj,kk)=feq + p_g(kk).*Hi ;%+ (1-omega_d).*(g(ii,jj,kk)-feq) ;
                 end
                 
            end
        end
    end
    for kk=1:5
        g(:,:,kk)=circshift(g(:,:,kk), [ex(kk),ey(kk),0]);
    end
    % bcs
    for ii=1:nx
        for jj=1:ny
            if(isfluid(ii,jj)==0)
                for kk=1:9
                    if(ii+ex(kk)>0 && jj+ey(kk)>0)
                        f(ii+ex(kk),jj+ey(kk),kk)=rho(ii,jj).*p(kk);
                    end
                end
                for kk=1:5
                    if(ii+ex(kk)>0 && jj+ey(kk)>0)
                        g(ii+ex(kk),jj+ey(kk),kk)=fi(ii,jj).*p_g(kk);
                    end
                end
                
            end
        end
    end
    fi(:,1)=fi(:,2);
    fi(:,ny)=fi(:,ny-1);
    fi(1,:)=fi(2,:);
    fi(nx,:)=fi(nx-1,:);
    %periodic x

    if(mod(tt,100)==0)      
       %imagesc(sqrt(u(2:end-1,2:end-1).^2+ v(2:end-1,2:end-1).^2)')
        %imagesc(sqrt(ffx.^2+ffy.^2)')
        %plot(u(80,2:end-1))
        imagesc(fi')
        %hold on
        %imagesc(mod_grad')
        %imagesc(rho)
        %contour(fi',[0.5 0.6])
        %imagesc(rho_phi(2:nx-1,2:ny-1)')
        axis xy
        axis equal
        colorbar
        %clim([0 max(max(sqrt(u(2:end-1,2:end-1).^2+ v(2:end-1,2:end-1).^2)'))])
        pause (0.01)
    end
    tt

end
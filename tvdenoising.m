% Rudin-Osher-Fatemi model
% INPUT u0 = noisy image, h=step size  lamda= parameter for denosing 
% OUTPUT: u(i,j) = denoised image 
% u0=u+n
% constraint of mean & std is overwhelmed by the property of method (no
% need to concern

function u=tvdenoising(u0,ITERATION,lamda)
u0=double(u0);
[M, N]=size(u0);
% u(x,y,0)=u0
u=zeros(M,N);
u(:,:)=u0;

% step size : h & #iteration N
% N_*h=1;
h=1/N;
% iterate N_ times
for Iter=1 : ITERATION
    %for every element of clear image u(i,j)
    %for i=2 : 1 : M-1
     % for j=2 : 1 :N-1
        UxM(2:M-1,2:N-1) = -u(1:M-2,2:N-1) + u(2:M-1,2:N-1);
        UxP(2:M-1,2:N-1) =  u(3:M,2:N-1)   - u(2:M-1,2:N-1);
        %UxM(i,j,Iter)=-(u(i-1,j,Iter)-u(i,j,Iter));
        %UxP(i,j,Iter)=(u(i+1,j,Iter)-u(i,j,Iter));
        %UxP0(i,j)=(u(i+1,j,1)-u(i,j,1));
        UyM(2:M-1,2:N-1) = -u(2:M-1,1:N-2) +u(2:M-1,2:N-1);
        UyP(2:M-1,2:N-1) =  u(2:M-1,3:N)   -u(2:M-1,2:N-1);
        %UyM(i,j,Iter)=-(u(i,j-1,Iter)-u(i,j,Iter));
        %UyP(i,j,Iter)=(u(i,j+1,Iter)-u(i,j,Iter));
        %UyP0(i,j)=(u(i,j+1,1)-u(i,j,1));
        m=@(a,b) min(abs(a),abs(b)).*(sign(a)+sign(b))/2;
    
        %lamda=-(h/2*std) * (sum(sum(hypot(UxP(i,j),UyP(i,j)) - (UxP0(i,j)*UxP(i,j)/(hypot(UxP(i,j),UyP(i,j))))  -(UyP0(i,j)*UyP(i,j))/(hypot(UxP(i,j),UyP(i,j))) )));
        Fco(2:M-1,2:N-1) = UxP(2:M-1,2:N-1)./(eps+(((UxP(2:M-1,2:N-1)).^2+m(UyP(2:M-1,2:N-1),UyM(2:M-1,2:N-1)).^2).^(1/2)));
        FcoxM(2:M-1,2:N-1)= - Fco(1:M-2,2:N-1) + Fco(2:M-1,2:N-1);
        %Fco(i,j,Iter)=UxP(i,j,Iter)./((UxP(i,j,Iter).^2 + UyP(i,j,Iter).^2).^(1/2));
        %FcoxM(i,j,Iter)=-(Fco(i-1,j,Iter)-Fco(i,j,Iter));
        Sco(2:M-1,2:N-1) = UyP(2:M-1,2:N-1)./(eps+(((UyP(2:M-1,2:N-1)).^2+m(UxP(2:M-1,2:N-1),UxM(2:M-1,2:N-1)).^2).^(1/2)));
        ScoyM(2:M-1,2:N-1)=-Sco(2:M-1,1:N-2) +Sco(2:M-1,2:N-1);
        %Sco(i,j,Iter)=UyP(i,j,Iter)./((UyP(i,j,Iter).^2 + UxP(i,j,Iter).^2).^(1/2));
        %ScoyM(i,j,Iter)=-(Sco(i,j-1,Iter)-Sco(i,j,Iter));
    
        derT=h*h/10;
        u(2:M-1,2:N-1)=u(2:M-1,2:N-1)+(derT/h)*(FcoxM(2:M-1,2:N-1)+ScoyM(2:M-1,2:N-1))-derT*lamda*(u(2:M-1,2:N-1)-u0(2:M-1,2:N-1));
        %u(i,j,Iter+1)=u(i,j,Iter) + ((derT/h)*(FcoxM(i,j,Iter)+ScoyM(i,j,Iter)- derT*lamda*(u(i,j,Iter)-u(i,j,1))));
        %difference(i,j,Iter)=u(i,j,Iter+1)-u(i,j,Iter);
        
%boundary condition : All derivative of edges are 0
%means that edges should be allocated by its neighborhood        

%Left & Right
    u(2:M-1,1)=u(2:M-1,2);
    u(2:M-1,N)=u(2:M-1,N);
%Up & down
    u(1,2:N-1)=u(2,2:N-1);
    u(M,2:N-1)=u(M-1,2:N-1); 

%Four pinpoint
    u(1,1)=u(2,2);
    u(1,N)=u(2,N-1); 
    u(M,1)=u(M-1,2);
    u(M,N)=u(M-1,N-1);
    if mod(Iter,ITERATION/5)==0
        figure, imshow(u(:,:),[])
    end
 
    energy(Iter) =sum(sum(hypot(UxP,UyP)))+ sum(sum(  (lamda/2)*(u-u0).^2  ));
end 
plot((1:ITERATION),(energy))
%result
figure,subplot(1,3,1), imshow(u0(:,:),[])
title('noisy image (u0)')

subplot(1,3,2) ,imshow(u(:,:),[])
title('clear image (u)')

subplot(1,3,3) ,imshow(u(:,:)-u0(:,:),[])
title('noise (n)')
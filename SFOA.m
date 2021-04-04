% SFOA optimization algorithm
function [jui,vsss]=SFOA
clc;
clear ;
tic



npop=pop();

maxit=round(itrc());

Fu=func();
[xmin,xmax,nvar,fobj]=Functions_info(Fu);

if size(xmin)==1
    
    Lb=xmin*ones(1,nvar);
    Ub=xmax*ones(1,nvar);
    
else
    
    Lb=xmin;
    Ub=xmax;
end
vmax=0.1*abs(Lb-Ub);


empty_sheep.position=[];
empty_sheep.velocity=[];
empty_sheep.cost=[];
empty_sheep.pbest=[];
empty_sheep.pbestcost=[];

sheep=repmat(empty_sheep,npop,1);
dxx=zeros(npop,nvar);
gbest=zeros(maxit,nvar);
gbestcost=zeros(maxit,1);

for it=1:maxit
    
    if it<=maxit*0.1
        gbestcost(1)=inf;
        if it>1
            gbest(it,:)=gbest(it-1,:);
            gbestcost(it)=gbestcost(it-1);
        end
        for i=1:npop
            kn=0;
            while(  kn==0   )
                sheep(i).velocity=zeros(1,nvar);
                x=Lb+(Ub-Lb).*rand(1,nvar);
                if it>1
                    for l=1:nvar
                        if rand>=0
                            sheep(i).position(l)=x(l);
                        end
                    end
                else
                    sheep(i).position=x;
                end
                
                if i<= npop*0.2
                    sheep(i).type=1;%goat
                else
                    sheep(i).type=0;%sheep
                end
                sheep(i).pbest=inf;
                sheep(i).pbestcost=inf;
                
                if ( Condition(  sheep(i).position,Fu)==1)
                    kn=1;
                end
                
            end
            sheep(i).cost=Cost(sheep(i).position);
            if sheep(i).cost<sheep(i).pbestcost
                sheep(i).pbest=sheep(i).position;
                sheep(i).pbestcost=sheep(i).cost;
                
                if sheep(i).pbestcost<gbestcost(it)
                    gbest(it,:)=sheep(i).pbest;
                    gbestcost(it)=sheep(i).pbestcost;
                end
            end
            sheep(i).position=   sheep(i).pbest;
            sheep(i).cost=    sheep(i).pbestcost;
        end
        
    else
        gbest(it,:)=gbest(it-1,:);
        gbestcost(it)=gbestcost(it-1);
        T=(1-it/maxit) ;
        
        for i=1:npop
            
            if mod(it,round(itrc/50))==0 || it>=maxit-5
                
                kn=0;
                while(  kn==0   )
                    if       dxx(i,:)==zeros(nvar,1)
                        if sheep(i).type>0.8
                            rg=(Ub-Lb)* (T) .* rand(1,nvar)*0.1;
                        else
                            rg=(Ub-Lb)* (T) .* rand(1,nvar)*0.001;
                        end
                        
                        dx=(2*rg)*rand-rg ;
                    else
                        dx=   dxx(i,:)*-1;
                        dxx(i,:)=zeros(nvar,1);
                    end
                    ps=sheep(i).position+dx;
                    
                    ps= min(max(ps,Lb),Ub);
                    
                    if ( Condition( ps,Fu)==1)
                        kn=1;
                    end
                end
                cs=Cost(ps);
                if  sheep(i).cost>cs
                    sheep(i).position=ps;
                    sheep(i).cost=cs;
                    
                    dxx(i,:)=zeros(nvar,1);
                    
                else
                    dxx(i,:)=dx;
                end
                
                
                
                
                if sheep(i).cost<sheep(i).pbestcost
                    sheep(i).pbest=sheep(i).position;
                    sheep(i).pbestcost=sheep(i).cost;
                    
                    if sheep(i).pbestcost<gbestcost(it)
                        gbest(it,:)=sheep(i).pbest;
                        gbestcost(it)=sheep(i).pbestcost;
                    end
                end
            else
                
                kn=0 ;
                while(  kn==0   )
                    
                    
                    if   sheep(i).type>=0.80
                        C=(3)*rand(1,nvar) ;
                        if T>0.7
                            
                            sheep(i).velocity=rand(1,nvar).*(gbest(it,:)-sheep(i).position)...
                                +(1-T).*2.*rand(1,nvar).*(sheep(i).pbest-sheep(i).position);
                        else
                            sheep(i).velocity=(1-T).*2.*rand(1,nvar).*gbest(it,:)-sheep(i).position;
                            
                        end
                        
                    else
                        C=(3)*rand ;
                        if T>0.3
                            randsheep=floor((npop)*rand);
                            randsheep=min(max(randsheep,1),npop);
                            sheep(i).velocity=C*rand(1,nvar).*(sheep(randsheep).pbest-sheep(i).position)...
                                +(1-T)*C*rand(1,nvar).*(gbest(it,:)-sheep(i).position)...
                                +C*rand(1,nvar).*(sheep(i).pbest-sheep(i).position);
                        else
                            sheep(i).velocity=C*rand(1,nvar).*(sheep(i).pbest-sheep(i).position)...
                                +C*(1-T).*(gbest(it,:)-sheep(i).position);
                            
                            
                        end
                        
                    end
                    
                    
                    
                    sheep(i).velocity=min(max(sheep(i).velocity,-vmax),vmax);
                    
                    for j=1:nvar
                        if rand<11
                            af=sheep(i).position;
                            af(j)=af(j)+sheep(i).velocity(j);
                            if Condition(af,Fu)==1
                                sheep(i).position(j)=sheep(i).position(j)+sheep(i).velocity(j);
                            end
                            
                            
                        end
                    end
                    sheep(i).position=min(max(sheep(i).position,Lb),Ub);
                    
                    sheep(i).position=(max(sheep(i).position,Lb));
                    if ( Condition(  sheep(i).position,Fu)==1)
                        kn=1;
                    end
                end
                sheep(i).cost=Cost(sheep(i).position);
                if sheep(i).cost<sheep(i).pbestcost
                    sheep(i).pbest=sheep(i).position;
                    sheep(i).pbestcost=sheep(i).cost;
                    
                    if sheep(i).pbestcost<gbestcost(it)
                        gbest(it,:)=sheep(i).pbest;
                        gbestcost(it)=sheep(i).pbestcost;
                    end
                end
                
            end
            
            
            
        end
        
    end
end
gh=zeros(1,itrc);
for i=1:round(itrc)
    gh(1,i)=  gbestcost(i,1);
end
jui=(gh);
vsss=gbest(itrc,:);
end

function G= Condition(x,Fu)
G=1;

end
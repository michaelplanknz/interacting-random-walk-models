using Plots, SpecialFunctions, Random, NLopt,Distributions, Interpolations, Dierckx, DifferentialEquations, LaTeXStrings
gr()
function Stochastic(LX,LY,t1,A0,PM,ρ)
Q=sum(A0)
T=Int(t1)
NumAgent=zeros(LX)
NumAgentT=zeros(LX)
    AAtemp=copy(A0)
    for i in 1:LX
        for j in 1:LY
           NumAgent[i]=NumAgent[i]+AAtemp[i,j]
       end
   end

    for kk in 1:T
        println(kk)
        count = 0
        while count  < Q
            II =rand(1:LX)
            JJ=rand(1:LY)
            α=rand(1)
            f=AAtemp[II,JJ]/maximum(AAtemp)
            if AAtemp[II,JJ] > 0.0 && II>1 && II<LX && JJ > 1 && JJ < LY &&  α[1] < f 
            count=count+1
            R =rand(1)
            S=rand(1)
                if R[1] > 0 && R[1]<=(1-ρ)/4 && S[1] <=PM
                AAtemp[II,JJ]=AAtemp[II,JJ]-1
                AAtemp[II-1,JJ]=AAtemp[II-1,JJ]+1
                elseif R[1] > (1-ρ)/4 && R[1]<=2/4  && S[1] <=PM
                AAtemp[II,JJ]=AAtemp[II,JJ]-1
                AAtemp[II+1,JJ]=AAtemp[II+1,JJ]+1
                elseif R[1] > 2/4 && R[1]<=3/4  && S[1] <=PM
                AAtemp[II,JJ]=AAtemp[II,JJ]-1
                AAtemp[II,JJ-1]=AAtemp[II,JJ-1]+1
                elseif R[1] > 3/4 && R[1]<=4/4  && S[1] <=PM
                AAtemp[II,JJ]=AAtemp[II,JJ]-1
                AAtemp[II,JJ+1]=AAtemp[II,JJ+1]+1
                end

                elseif AAtemp[II,JJ] > 0.0 && II>1 && II<LX && JJ == 1 &&  α[1] < f 
                    count=count+1
                    R =rand(1)
                    S=rand(1)
                        if R[1] > 0 && R[1]<=(1-ρ)/4  && S[1] <=PM
                        AAtemp[II,JJ]=AAtemp[II,JJ]-1
                        AAtemp[II-1,JJ]=AAtemp[II-1,JJ]+1
                        elseif R[1] > (1-ρ)/4  && R[1]<=2/4  && S[1] <=PM
                        AAtemp[II,JJ]=AAtemp[II,JJ]-1
                        AAtemp[II+1,JJ]=AAtemp[II+1,JJ]+1
                        elseif R[1] > 3/4 && R[1]<=4/4 && S[1] <=PM
                        AAtemp[II,JJ]=AAtemp[II,JJ]-1
                        AAtemp[II,JJ+1]=AAtemp[II,JJ+1]+1
                        end
        

                    elseif AAtemp[II,JJ] > 0.0 && II>1 && II<LX &&  JJ == LY &&  α[1] < f 
                            count=count+1
                            R =rand(1)
                            S=rand(1)
                                if R[1] > 0 && R[1]<=(1-ρ)/4   && S[1] <=PM
                                AAtemp[II,JJ]=AAtemp[II,JJ]-1
                                AAtemp[II-1,JJ]=AAtemp[II-1,JJ]+1
                                elseif R[1] > (1-ρ)/4  && R[1]<=2/4  && S[1] <=PM
                                AAtemp[II,JJ]=AAtemp[II,JJ]-1
                                AAtemp[II+1,JJ]=AAtemp[II+1,JJ]+1
                                elseif R[1] > 2/4 && R[1]<=3/4 && S[1] <=PM
                                AAtemp[II,JJ]=AAtemp[II,JJ]-1
                                AAtemp[II,JJ-1]=AAtemp[II,JJ-1]+1
                                 end
            end
        end
    end


for i in 1:LX
     for j in 1:LY
        NumAgentT[i]=NumAgentT[i]+AAtemp[i,j]
    end
end

return NumAgent./LY, NumAgentT./LY
end



function diff!(du,u,p,t)
    dx,N,D,v=p 
    for i in 2:N-1
    du[i]=D*(u[i-1]-2*u[i]+u[i+1])/dx^2-v*(u[i+1]-u[i-1])/(2*dx)
    end
    du[1]=D*(u[2]-u[1])/dx^2  
    du[N]=D*(u[N-1]-u[N])/dx^2
    end
      
        
    
    
    function pdesolver(LX,dx,N,T,h,U0,D,v)
    p=(dx,N,D,v)
    ic=zeros(N)
    for i in 1:N
        if abs(-LX/2+(i-1)*dx) <= h
        ic[i]=U0
        end
    end    
    tspan=(0.0,T)
    prob=ODEProblem(diff!,ic,tspan,p)
   sol=solve(prob,saveat=T);
    return sol;
    end








LX=200
LY=100
P=1.0
ρ=0.5
D=P/4
v=P*ρ/2
T=200
U0=1.0
h=10
dx=0.5;
N=Int(round(LX/dx))+1;

A0=zeros(LX,LY)
xxloc=zeros(LX)
yyloc=zeros(LY)
for i in 1:LX
    xxloc[i]=-LX/2+(i-1)
    for j in 1:LY
    R=rand()
    yyloc[j]=0+(j-1)
        if abs(xxloc[i]) <= h && R <= U0
        A0[i,j]=1.0
        end
    end
end

(density0,densityT)=Stochastic(LX,LY,T,A0,P,ρ);
numsol=pdesolver(LX,dx,N,T,h,U0,D,v);
#p1=plot(xxloc,density0,label=false,lw=4,lc=:red)
p1=plot(-LX/2:dx:LX/2,numsol[:,1],label=false,lw=2,lc=:black)
p1=plot!(-LX/2:dx:LX/2,numsol[:,2],label=false,lw=2,lc=:black)
p1=plot!(xxloc,densityT,label=false,lw=2,lc=:red)
q1=plot!(xlims=(-100,100),xticks=([-100,-50,0,50,100],[L"-100",L"-50",L"0",L"50",L"100"]))
q1=plot!(ylims=(0,1.00),yticks=([0,0.25,0.5,0.75,1.0],[L"0.00", L"0.25", L"0.50", L"0.75", L"1.00"]))
q1=plot!(xguidefontsize=12, yguidefontsize=12,xtickfontsize=12, ytickfontsize=12)
q1=plot!(xlabel=L"x",ylabel=L"u(x,t)")
q1=plot!(title="(b) Biased noninteracting")
savefig(q1,"Figure.pdf")





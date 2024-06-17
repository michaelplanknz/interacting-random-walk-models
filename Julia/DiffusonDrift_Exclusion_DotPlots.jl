using Plots, SpecialFunctions, Random, NLopt,Distributions, Interpolations, Dierckx, LaTeXStrings
gr()
function Stochastic(Q,LX,LY,t1,A0,PM,ρ,tk1,tk2,tk3) 
    pos0=zeros(Q,2)
    post=zeros(Q,2)
    T=Int(t1)
    AA=zeros(LX,LY)
    NumAgent=zeros(LX)
    NumAgentT=zeros(LX)
        AAtemp=copy(A0)
        
        agent=0
        for i in 1:LX
            for j in 1:LY
                if A0[i,j] > 0.0
                agent = agent+1
                pos0[agent,1]=i
                pos0[agent,2]=j
                end 
            end
          end

    
        for i in 1:LX
            for j in 1:LY
               NumAgent[i]=NumAgent[i]+AAtemp[i,j]
           end
       end
    
        for kk in 1:T
            count = 0
            while count  < Q
                II =rand(1:LX)
                JJ=rand(1:LY)
                if AAtemp[II,JJ] > 0.0 && II>1 && II<LX && JJ > 1 && JJ < LY 
                count=count+1
                R =rand(1)
                S=rand(1)
                    if AAtemp[II-1,JJ] == 0.0 && R[1] > 0 && R[1]<=(1-ρ)/4 && S[1] <=PM
                    AAtemp[II-1,JJ]=AAtemp[II,JJ]
                    AAtemp[II,JJ]=0.0
                    elseif AAtemp[II+1,JJ] == 0.0 && R[1] > (1-ρ)/4 && R[1]<=2/4  && S[1] <=PM
                    AAtemp[II+1,JJ]=AAtemp[II,JJ]
                    AAtemp[II,JJ]=0.0
                    elseif AAtemp[II,JJ-1] == 0.0 && R[1] > 2/4 && R[1]<=3/4  && S[1] <=PM
                    AAtemp[II,JJ-1]=AAtemp[II,JJ]
                    AAtemp[II,JJ]=0.0
                    elseif AAtemp[II,JJ+1] == 0.0 && R[1] > 3/4 && R[1]<=4/4  && S[1] <=PM
                    AAtemp[II,JJ+1]=AAtemp[II,JJ]
                    AAtemp[II,JJ]=0.0
                    end
    
                    elseif AAtemp[II,JJ] > 0.0 && II>1 && II<LX && JJ == 1 
                        count=count+1
                        R =rand(1)
                        S=rand(1)
                        if AAtemp[II-1,JJ] == 0.0 && R[1] > 0 && R[1]<=(1-ρ)/4 && S[1] <=PM
                            AAtemp[II-1,JJ]=AAtemp[II,JJ]
                            AAtemp[II,JJ]=0.0
                            elseif AAtemp[II+1,JJ] == 0.0 && R[1] > (1-ρ)/4 && R[1]<=2/4  && S[1] <=PM
                            AAtemp[II+1,JJ]=AAtemp[II,JJ]
                            AAtemp[II,JJ]=0.0
                            elseif AAtemp[II,JJ+1] == 0.0 && R[1] > 3/4 && R[1]<=4/4  && S[1] <=PM
                            AAtemp[II,JJ+1]=AAtemp[II,JJ]
                            AAtemp[II,JJ]=0.0
                            end
            
    
                        elseif AAtemp[II,JJ] > 0.0 && II>1 && II<LX &&  JJ == LY 
                                count=count+1
                                R =rand(1)
                                S=rand(1)
                                if AAtemp[II-1,JJ] == 0.0 && R[1] > 0 && R[1]<=(1-ρ)/4 && S[1] <=PM
                                    AAtemp[II-1,JJ]=AAtemp[II,JJ]
                                    AAtemp[II,JJ]=0.0
                                    elseif AAtemp[II+1,JJ] == 0.0 && R[1] > (1-ρ)/4 && R[1]<=2/4  && S[1] <=PM
                                    AAtemp[II+1,JJ]=AAtemp[II,JJ]
                                    AAtemp[II,JJ]=0.0
                                    elseif AAtemp[II,JJ-1] == 0.0 && R[1] > 2/4 && R[1]<=3/4  && S[1] <=PM
                                    AAtemp[II,JJ-1]=AAtemp[II,JJ]
                                    AAtemp[II,JJ]=0.0
                                    end
                end
            end
        for i in 1:LX  
            for j in 1:LY
                if AAtemp[i,j] == 10
                    tk1[1,kk+1]=(i-1)
                    tk1[2,kk+1]=(j-1)
                elseif AAtemp[i,j]==20
                    tk2[1,kk+1]=(i-1)
                    tk2[2,kk+1]=(j-1)
                elseif AAtemp[i,j]==30
                    tk3[1,kk+1]=(i-1)
                    tk3[2,kk+1]=(j-1)    
                end    
                
            end
        end







        end
    AA=AA+AAtemp
    
    for i in 1:LX
         for j in 1:LY
            NumAgentT[i]=NumAgentT[i]+AA[i,j]
        end
    end
    

    agent=0
    for i in 1:LX
        for j in 1:LY
            if AA[i,j] == 1.0
            agent = agent+1
            post[agent,1]=i
            post[agent,2]=j
            end 
        end
      end


    return NumAgent./LY, NumAgentT./LY,pos0,post,tk1,tk2,tk3
    end




LX=200
LY=50
PM=1.0
ρ=0.5
T=200.0
U0=1.0
h=10
Q=0
A0=zeros(LX,LY)
xxloc=zeros(LX)
yyloc=zeros(LY)
for i in 1:LX
    xxloc[i]=-LX/2+(i-1)
    for j in 1:LY
    R=rand()
    yyloc[j]=0+(j-1)
        if abs(xxloc[i]) <= h && R <= U0
            Q=Q+1
            A0[i,j]=1.0
        end
    end
end
A0[111,26]=10
A0[91,26]=20
A0[101,26]=30
tk1=zeros(2,Int(T)+1)
tk1[1,1]=(111-1)
tk1[2,1]=(26-1)
tk2=zeros(2,Int(T)+1)
tk2[1,1]=(91-1)
tk2[2,1]=(26-1)
tk3[1,1]=(101-1)
tk3[2,1]=(26-1)


(data0,data,pos0,post,tk1,tk2,tk3)=Stochastic(Q,LX,LY,T,A0,PM,ρ,tk1,tk2,tk3)





q1=scatter(pos0[:,1],pos0[:,2],markersize=1.5,markershape=:circle, markercolor=:red,xlims=(0,LX),ylims=(0,LY+0.1),legend=false,msw=0,aspect_ratio=:equal)
q1=plot!(xticks = ([1,50,100,150,200],  [L"-100", L"-50", L"0", L"50", L"100"]))
q1=plot!(yticks = ([0,25,50],[L"0", L"25", L"50"]))
q1=plot!(xlabel=L"x",ylabel=L"y")
q1=plot!(xguidefontsize=12, yguidefontsize=12,xtickfontsize=12, ytickfontsize=12)
q1=plot!(title="(d) Initial condition")
q2=scatter(post[:,1],post[:,2],markersize=1.5,markershape=:circle, markercolor=:red,xlims=(0,LX),ylims=(0,LY+0.1),legend=false,msw=0,aspect_ratio=:equal)
q2=plot!(xticks = ([1,50,100,150,200],  [L"-100", L"-50", L"0", L"50", L"100"]))
q2=plot!(yticks = ([0,25,50],[L"0", L"25", L"50"]))
q2=plot!(xlabel=L"x",ylabel=L"y")
q2=plot!(xguidefontsize=12, yguidefontsize=12,xtickfontsize=12, ytickfontsize=12)
q2=plot!(tk1[1,:],tk1[2,:],lw=3,lc=:green)
q2=plot!(tk2[1,:],tk2[2,:],lw=3,lc=:blue)
q2=plot!(tk3[1,:],tk3[2,:],lw=3,lc=:turquoise1)
q2=plot!(title="(f) Biased interacting")
q3=plot(q1,q2,layout=(2,1))
display(q3)
savefig(q3,"Figure.pdf")
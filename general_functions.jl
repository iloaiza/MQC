function super_zeros(n,dims,le_type="Float64")
    zero_string="Z=[zeros($le_type,$n) for i1 in 1:$(dims[1])"
    for i in 2:length(dims)
        zero_string=zero_string*", i$i in 1:$(dims[2])"
    end
    zero_string=zero_string*"]"
    eval(Meta.parse(zero_string))

    return Z
end

function NxN(BASE,ex,N=length(BASE),le_type=Float64)
    #build (res)x(res)x...x(res) (N times) matrix with N dimensional vectors in each entry
    #BASE is an array of linspaces over which the array is built
    #ex is the number of extra dimensions to add over N
    #le_type admits string entry, change to "Complex" for storing complex stuff
    #or maybe even touples ;)
    global A_NxN=super_zeros(N+ex,length.(BASE),le_type)
    global B_NxN=BASE
    #@show size(A)
    forstring="for i1 in 1:$(length(BASE[1]))"
    for i in 2:N
        forstring=forstring*", i$i in 1:$(length(BASE[i]))"
    end
    forstring=forstring*"    A_NxN[i1"
    for i in 2:N
        forstring=forstring*",i$i"
    end
    forstring=forstring*"]=[B_NxN[1][i1]"
    for i in 2:N
        forstring=forstring*",B_NxN[$i][i$i]"
    end
    for i in 1:ex
        forstring=forstring*",0"
    end
    forstring=forstring*"];  end"
    eval(Meta.parse(forstring))

    return A_NxN
end


function my_histogram(X,n)
    Xrange=range(minimum(X),stop=maximum(X),length=n+1)
    plotting_X=range(Xrange[1]+step/2,stop=Xrange[end]-step/2,length=n)
    H=zeros(n)

    for i in 1:n-1
        H[i]=count(x->(Xrange[i]<=x<Xrange[i+1]),X)
    end
    H[n]=count(x->(Xrange[n]<=x<=Xrange[n+1]),X)

    return plotting_X,H
end

function my_histogram(X,xmin,xmax,n)
    Xrange=range(xmin,stop=xmax,length=n+1)
    step=Xrange[2]-Xrange[1]
    plotting_X=range(Xrange[1]+step/2,stop=Xrange[end]-step/2,length=n)
    H=zeros(n)

    for i in 1:n-1
        H[i]=count(x->(Xrange[i]<=x<Xrange[i+1]),X)
    end
    H[n]=count(x->(Xrange[n]<=x<=Xrange[n+1]),X)

    return plotting_X,H
end

function super_histo(R,n)
    xmin=minimum(R)
    xmax=maximum(R)
    flags=length(R[:,1])
    plotting_R=zeros(n,flags)
    histo_R=zeros(n,flags)
    for i in 1:flags
        plotting_R[:,i],histo_R[:,i]=my_histogram(R[i,:],xmin,xmax,n)
    end

    return plotting_R,histo_R
end

function super_histo(R,xmin,xmax,n)
    flags=length(R[:,1])
    plotting_R=zeros(n,flags)
    histo_R=zeros(n,flags)
    for i in 1:flags
        plotting_R[:,i],histo_R[:,i]=my_histogram(R[i,:],xmin,xmax,n)
    end

    return plotting_R,histo_R
end

function integral(f,a,b,dx=1e-5) #Riemann integral
    X=collect(a:dx:b)
    SUM=0
    for x in X
        SUM+=f(x)
    end
    SUM=SUM*dx

    return SUM
end

function eye(dims,eltype=Float64)
    M=zeros(eltype,dims,dims)
    for n in 1:dims
        M[n,n]=1
    end

    return M
end

function norm(v)
    return sqrt(sum(abs2.(v)))
end

function norm2(v)
    return sum(abs2.(v))
end

function KEP = Approximate_keppler_elements(data,GM)
    r_icrs=data(:,4:6);
    t=data(:,3);

    rt_icrs=(r_icrs(2,:)-r_icrs(1,:))/(t(2)-t(1));
    r_icrs=r_icrs(1,:);
    t=t(1);

    h=cross(r_icrs,rt_icrs);
    W=atan2(h(1),-h(2));

    h_=norm(h);
    i=acos(h(3)/h_);

    r=norm(r_icrs);
    A=(dot(r_icrs,rt_icrs)*h_)/(GM*r);
    B=(h_^2)/(GM*r)-1;
    e=sqrt(A^2+B^2);

    a=(h_^2)/(GM*(1-e^2));

    v=atan2(A,B);
    E=2*atan(sqrt((1-e)/(1+e))*tan(v/2));
    M=E-e*sin(E);
    n=sqrt(GM/(a^3));
    t0=t-(M/n);
    
    rxy_c =   [1       0        0
               0   cos(i)   sin(i)
               0  -sin(i)   cos(i)]*[cos(W)   sin(W)     0
                                    -sin(W)   cos(W)     0
                                         0        0      1]*r_icrs';
    u=atan2(rxy_c(2),rxy_c(1));
    w=u-v;

    
    KEP=[a,e,t0,w,i,W];
end